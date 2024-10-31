def open_file(file_path):
    import gzip
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file

def confidence_interval(line, lower =0.01, upper = 0.98):
    """
    retain variants for which the confidence interval of variant allele frequency (VAF) overlaps the range 1% to 98%
    """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    GT,DP,HF,CILOW,CIUP,SDP = sample.strip().split(":")
    if lower <= float(CILOW) <= float(CIUP) <= upper:
        return True
    else:
        return False
    
def rm_indel_mul_alt(line):
    """
    remove all indels, as detecting heteroplasmic indels is unreliable using short-read sequencing data; 
    remove variants at sites with multiple alternate alleles
    """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    if len(ref) == len(alt) == 1:
        return True
    else:
        return False

def rm_region(line, region_list=[[66,71], [300,316], [513,525], [3106,3107], [12418,12425],[16182,16194]]):
    """
    remove variants falling within specific regions associated with misalignment errors related to homopolymeric tracts (np 66-71, 300-316, 513-525, 3106-3107, 12418-12425 and 16182-16194)
    """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    poses = []
    for list in region_list:
        for i in range(list[0], list[1]+1):
            poses.append(i)
    if int(pos) in poses:
        return False
    else:
        return True
    
def rm_vars_lt2reads(line):
    """
    6) remove variants less than 2 reads on each strand with the minor allele; 
    """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    GT,DP,HF,CILOW,CIUP,SDP = sample.strip().split(":")
    ws_count, cs_count = SDP.strip().split(';')
    if float(ws_count) < 2 or float(cs_count) <2:
        return False
    else:
        return True

def rm_depthlt200_HF(line):
    """
    7) remove heteroplasmic variants with depth < 200x and low level heteroplasmies (HFs<5%) with depth < 500x,
    ## unless a Bonferroni-corrected p-value < 10-5 was given by deepSNV 
    """
    chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
    GT,DP,HF,CILOW,CIUP,SDP = sample.strip().split(":")
    if float(DP)<200:
        return False
    elif float(HF) < 0.05 and float(DP) < 500:
        return False
    else:
        return True
        

def filter_vcf(vcf_file, out_file_just_add_LowQual, out_file_filtered):
    with open(out_file_just_add_LowQual, 'w') as out_file_add , open(out_file_filtered, 'w') as out_file_rm:
        vcff = open_file(vcf_file)
        out_file_add.write("##fileformat=VCFv4.2\n##FILTER=<ID=LowQual,Description='Low Quality'>\n")
        out_file_rm.write("##fileformat=VCFv4.2\n##FILTER=<ID=LowQual,Description='Low Quality'>\n")
        for line in vcff:
            if line.startswith('#'):
                out_file_add.write(line)
                out_file_rm.write(line)
            else:
                chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
                if rm_indel_mul_alt(line):
                    if (confidence_interval(line) and rm_region(line) 
                        and rm_vars_lt2reads(line) and rm_depthlt200_HF(line)):
                        out_file_rm.write(line)
                        out_file_add.write(line)
                    else:
                        filter = 'LowQual'
                        print(chrom, pos, id, ref, alt, qual, filter, info, format, sample, sep='\t', file=out_file_add)
                else:
                    filter = 'LowQual'
                    print(chrom, pos, id, ref, alt, qual, filter, info, format, sample, sep='\t', file=out_file_add)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='filter vcf file')
    parser.add_argument('-i', '--input', required=True, help='input vcf file')
    parser.add_argument('-o', '--output', required=True, help='output filtered vcf file')
    parser.add_argument('-l', '--lowqual', required=True, help='output vcf file with LowQual information')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    filter_vcf(args.input, args.lowqual, args.output)

