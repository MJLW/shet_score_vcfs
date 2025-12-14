use std::borrow::Cow;
use std::io::BufWriter;
use std::iter::zip;
use std::{error::Error, fs::File, io::BufReader, path::Path};

use clap::Parser;
use log::{error, warn};
use noodles::vcf::record::Info;
use noodles::vcf::variant::record::info::field::{Value, value::Array};
use noodles::vcf::variant::record::samples::Series;
use noodles::vcf::variant::record::samples::keys::key;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
use noodles::{bgzf, vcf};
use regex::Regex;
use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(name = "", version = env!("CARGO_PKG_VERSION"), about = "", long_about = "")]
struct Args {
    #[arg(short, long)]
    shet_scores: String,

    #[arg(short, long)]
    vcfs_file: String,

    #[arg(short, long)]
    out: String,
}

#[derive(Deserialize)]
struct Row {
    #[serde(rename = "ID")]
    id: String,

    #[serde(rename = "VCF")]
    vcf_path: String,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct GeneSHET {
    #[serde(rename = "Gene")]
    gene: String,

    #[serde(rename = "U")]
    u: f64,

    #[serde(rename = "n_NFE")]
    n_nfe: f64,

    #[serde(rename = "k_NFE")]
    k_nfe: f64,

    s_het_det: f64,
    low_det: f64,
    up_det: f64,
    s_het_drift: f64,
    low_drift: f64,
    up_drift: f64,
}

#[derive(Serialize)]
struct Output {
    id: String,
    paternal_haplotype: f64,
    maternal_haplotype: f64,
}

pub trait InfoValue {
    fn as_f32(&self) -> Result<f32, Box<dyn Error>>;
    fn as_string(&self) -> Result<Cow<'_, str>, Box<dyn Error>>;
    fn as_i32_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<i32>, std::io::Error>> + 'a>>;

    fn as_f32_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<f32>, std::io::Error>> + 'a>>;
    fn as_string_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<String>, std::io::Error>> + 'a>>;
}

impl InfoValue for Value<'_> {
    fn as_f32(&self) -> Result<f32, Box<dyn Error>> {
        match self {
            Value::Float(float_value) => {
                println!("{}", float_value.to_owned());
                return Ok(float_value.to_owned());
            }
            _ => {
                Err("Could not parse tag to value. It does not match the expected encoding or is unexpectedly missing.".into())
            }
        }
    }

    fn as_string(&self) -> Result<Cow<'_, str>, Box<dyn Error>> {
        match self {
            Value::String(string_value) => return Ok(string_value.to_owned()),
            _ => {
                Err("Could not parse tag to value. It does not match the expected encoding or is unexpectedly missing.".into())
            }
        }
    }

    fn as_i32_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<i32>, std::io::Error>> + 'a>> {
        match self {
            Value::Array(array_value) => match array_value {
                Array::Integer(integer_values) => Some(Box::new(integer_values.iter())),
                _ => None,
            },
            _ => None,
        }
    }

    fn as_f32_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<f32>, std::io::Error>> + 'a>> {
        match self {
            Value::Array(array_value) => match array_value {
                Array::Float(float_values) => Some(Box::new(float_values.iter())),
                _ => None,
            },
            _ => None,
        }
    }

    // Could not figure out how to return a Cow<'a, str> object, the lifetimes are screwy, so
    // instead it returns a String
    fn as_string_array<'a>(
        &'a self,
    ) -> Option<Box<dyn Iterator<Item = Result<Option<String>, std::io::Error>> + 'a>> {
        match self {
            Value::Array(array_value) => match array_value {
                Array::String(string_values) => {
                    Some(Box::new(string_values.iter().map(|s| match s {
                        Ok(Some(cow)) => Ok(Some(cow.to_string())),
                        Ok(None) => Ok(None),
                        Err(e) => Err(e),
                    })))
                }
                _ => None,
            },
            _ => None,
        }
    }
}

fn expect_value<T>(value: Option<T>, tag: &str) -> T {
    value.expect(&format!(
        "Missing value in INFO/{}, these are required for this program.",
        tag
    ))
}

fn read_input(path: &str) -> Result<Vec<Row>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(Path::new(path))?));

    let rows: Vec<Row> = reader
        .deserialize()
        .into_iter()
        .filter_map(|row| match row {
            Ok(r) => Some(r),
            Err(e) => {
                error!("Skipped row due to failed parsing: {}", e);
                None
            }
        })
        .collect();

    Ok(rows)
}

fn read_shet(path: &str) -> Result<Vec<GeneSHET>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(Path::new(path))?));

    let rows: Vec<GeneSHET> = reader
        .deserialize()
        .into_iter()
        .filter_map(|row| match row {
            Ok(r) => Some(r),
            Err(e) => {
                error!("Skipped row due to failed parsing: {}", e);
                None
            }
        })
        .collect();

    Ok(rows)
}

fn get_tag_value_from_info<'r>(
    info: &'r Info,
    header: &'r vcf::Header,
    tag: &str,
) -> Option<Value<'r>> {
    let tag_value = info.get(header, tag);

    return match tag_value {
        Some(Ok(Some(tag))) => Some(tag),
        _ => None,
    };
}

fn get_gt(record: &vcf::Record, header: &vcf::Header) -> Result<Genotype, Box<dyn Error>> {
    let samples = record.samples();
    let gt_series = samples.select(key::GENOTYPE).expect("Missing GT field");

    let gt_value = gt_series
        .get(header, 0)
        .expect("No GT field.")
        .expect("No GT value.");

    let Ok(vcf::variant::record::samples::series::Value::Genotype(genotype)) = gt_value else {
        panic!("Failed to cast genotype series element, underlying data is invalid.");
    };

    Ok(genotype.as_ref().try_into().unwrap())
}

fn write_output(path: &str, output_row: Vec<Output>) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(BufWriter::new(File::create(Path::new(path))?));

    // Write out header
    writer.write_record(&["ID", "PATERNAL_HAPLOTYPE", "MATERNAL_HAPLOTYPE"])?;

    for row in output_row {
        writer.serialize((row.id, row.paternal_haplotype, row.maternal_haplotype))?;
    }

    writer.flush()?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    const AF_TAG: &str = "MAX_AF";
    const AF_MIN_VALUE: f32 = 0.001;
    const VEP_CONSEQUENCE_TAG: &str = "vepConsequence";
    const VEP_BIOTYPE: &str = "vepBIOTYPE";
    const RESCUE_TAG: &str = "RESCUE";
    const PTV_TAG: &str = "PTV";
    const VEP_GENE_SYMBOL: &str = "vepSYMBOL";

    let ptv_csq_re =
        Regex::new("stop_gained|frameshift_variant|splice_acceptor_variant|splice_donor_variant")?;

    pretty_env_logger::init();
    let args = Args::parse();

    let shet_scores = read_shet(&args.shet_scores)?;

    let vcf_files = read_input(&args.vcfs_file)?;

    let mut output_rows: Vec<Output> = vec![];
    for row in vcf_files {
        let mut vcf = File::open(row.vcf_path)
            .map(bgzf::io::Reader::new)
            .map(vcf::io::Reader::new)?;
        let vcf_header = vcf.read_header()?;

        let mut paternal: f64 = 0.0;
        let mut maternal: f64 = 0.0;
        for record_result in vcf.records() {
            let record = record_result?;
            let info = record.info();

            // Check only rare variants
            // Gets the first value out, multiallelics are split
            let af_info_value = get_tag_value_from_info(&info, &vcf_header, AF_TAG).unwrap();

            // The expect_value function should probably be handled as .expect(custom error given tag) instead
            let mut af_info_array = expect_value(af_info_value.as_f32_array(), AF_TAG);
            let af = expect_value(expect_value(af_info_array.next(), AF_TAG)?, AF_TAG);

            if af > AF_MIN_VALUE {
                continue;
            }

            // Check if record is a PTV in VEP
            let vep_consequence_value =
                get_tag_value_from_info(&info, &vcf_header, VEP_CONSEQUENCE_TAG).unwrap();
            let vep_consequence_array =
                expect_value(vep_consequence_value.as_string_array(), VEP_CONSEQUENCE_TAG);

            let is_ptv: Vec<bool> = vep_consequence_array
                .into_iter()
                .map(|consequence| match consequence {
                    Ok(Some(csq)) => csq.split('&').any(|c| ptv_csq_re.is_match(c)),
                    _ => false,
                })
                .collect();

            // Lots of variants will not be a PTV is VEP, so we can escape those here
            let any_ptv = is_ptv.iter().any(|b| *b);
            if !any_ptv {
                continue;
            }

            // Make sure the variant is protein coding, otherwise there was never anything ran on them
            let biotype_value = get_tag_value_from_info(&info, &vcf_header, VEP_BIOTYPE).unwrap();
            let biotype_array = expect_value(biotype_value.as_string_array(), VEP_BIOTYPE);

            let is_protein_coding: Vec<bool> = biotype_array
                .into_iter()
                .map(|biotype| match biotype {
                    Ok(Some(biotype)) => biotype == "protein_coding",
                    _ => false,
                })
                .collect();

            let any_protein_coding = is_protein_coding
                .iter()
                .zip(is_ptv.clone())
                .any(|(protein_coding, vep_ptv)| *protein_coding && vep_ptv);
            if !any_protein_coding {
                continue;
            }

            // Figure out which genotypes rescue should be checked for
            let gt: Genotype = get_gt(&record, &vcf_header)?;
            let paternal_gt: &Allele = gt
                .as_ref()
                .get(0)
                .expect("Failed to find paternal haplotype genotype.");
            let maternal_gt: &Allele = gt
                .as_ref()
                .get(1)
                .expect("Failed to find maternal haplotype genotype.");

            let is_phased: bool = paternal_gt.phasing() == Phasing::Phased
                || maternal_gt.phasing() == Phasing::Phased;

            let paternal_call = paternal_gt
                .position()
                .expect("Could not find paternal haplotype call");
            let maternal_call = maternal_gt
                .position()
                .expect("Could not find maternal haplotype call");

            // If unphased, then heterozygous variants were applied to both haplotypes in RescueRanger
            let possible_paternal_call = (paternal_call > 0) || (!is_phased && (maternal_call > 0));
            let possible_maternal_call = (maternal_call > 0) || (!is_phased && (paternal_call > 0));

            // println!(
            //     "{}\t{}\t{}",
            //     row.id,
            //     record.reference_sequence_name(),
            //     record.variant_start().unwrap()?
            // );

            // Check RR PTVs
            let rescue_value = get_tag_value_from_info(&info, &vcf_header, RESCUE_TAG);
            if rescue_value.is_none() {
                continue;
            }
            let rescue_value_unwrapped = rescue_value.unwrap();
            let rescue_array_wrapped = rescue_value_unwrapped.as_string_array();

            let rescue_array = expect_value(rescue_array_wrapped, RESCUE_TAG);
            let rescues: Vec<_> = rescue_array
                .into_iter()
                .map(|rescue| match rescue {
                    Ok(Some(r)) => {
                        let rescue_haplotypes = r.split_at(1);
                        let paternal = rescue_haplotypes.0.parse::<i64>().unwrap_or(0);
                        let maternal = rescue_haplotypes.1.parse::<i64>().unwrap_or(0);
                        Some((paternal, maternal))
                    }
                    _ => None,
                })
                .collect();

            let ptv_value = get_tag_value_from_info(&info, &vcf_header, PTV_TAG).unwrap();
            let ptv_array = expect_value(ptv_value.as_string_array(), PTV_TAG);
            let ptvs: Vec<_> = ptv_array
                .into_iter()
                .map(|ptv| match ptv {
                    Ok(Some(r)) => {
                        let ptv_haplotypes = r.split_at(1);
                        let paternal = ptv_haplotypes.0.parse::<i64>().unwrap_or(0);
                        let maternal = ptv_haplotypes.1.parse::<i64>().unwrap_or(0);
                        Some((paternal, maternal))
                    }
                    _ => None,
                })
                .collect();

            // Find genes corresponding to PTV notations
            let gene_value = get_tag_value_from_info(&info, &vcf_header, VEP_GENE_SYMBOL).unwrap();
            let gene_array = expect_value(gene_value.as_string_array(), VEP_GENE_SYMBOL);

            let genes: Vec<_> = gene_array
                .into_iter()
                .map(|gene| match gene {
                    Ok(Some(gene)) => Some(gene),
                    _ => None,
                })
                .collect();

            // Get gene symbols for unrescued variants
            let mut affected_genes: Vec<String> = vec![];
            let mut affected_haplotypes: Vec<(bool, bool)> = vec![];
            for i in 0..genes.len() {
                let is_ptv = *is_ptv.get(i).unwrap() && *is_protein_coding.get(i).unwrap();

                if !is_ptv {
                    continue;
                }

                let mut paternal_rescued = false;
                let mut maternal_rescued = false;
                if let Some(rescue) = rescues.get(i) {
                    if let Some(rescue_tuple) = rescue {
                        paternal_rescued = rescue_tuple.0 > 0 && possible_paternal_call;
                        maternal_rescued = rescue_tuple.1 > 0 && possible_maternal_call;
                    }
                }

                let mut paternal_ptv = false;
                let mut maternal_ptv = false;
                if let Some(ptv) = ptvs.get(i) {
                    if let Some(ptv_tuple) = ptv {
                        paternal_ptv = ptv_tuple.0 > 0 && possible_paternal_call;
                        maternal_ptv = ptv_tuple.1 > 0 && possible_maternal_call;
                    }
                }

                let paternal_affected = !paternal_rescued || paternal_ptv;
                let maternal_affected = !maternal_rescued || maternal_ptv;

                if !paternal_affected && !maternal_affected {
                    continue;
                }

                if let Some(gene) = genes.get(i) {
                    let gene_name = gene.clone().unwrap();
                    if affected_genes.contains(&gene_name) {
                        continue;
                    }

                    affected_genes.push(gene_name);
                    affected_haplotypes.push((paternal_affected, maternal_affected));
                };
            }

            for (gene, haplotypes) in zip(affected_genes, affected_haplotypes) {
                let shet_search = shet_scores.iter().find(|shet| shet.gene == gene);

                if shet_search.is_none() {
                    warn!("Could not find SHET score for gene {}.", gene);
                    continue;
                }

                let matched_shet = shet_search.unwrap();

                if haplotypes.0 {
                    paternal += matched_shet.s_het_drift;
                }

                if haplotypes.1 {
                    maternal += matched_shet.s_het_drift;
                }
            }
        }

        let output_row: Output = Output {
            id: row.id.clone(),
            paternal_haplotype: paternal,
            maternal_haplotype: maternal,
        };

        output_rows.push(output_row);
    }

    write_output(&args.out, output_rows)?;

    Ok(())
}
