from collections import defaultdict
import gzip
import argparse

def genome2str(genome):
    """
    Convert a genome file to a string.
    """
    str = ''
    with open(genome, 'r') as f:
        for line in f:
            if line.startswith('>'):
                pass
            else:
                str += line.strip()
    return str

def coverage2dict(coverage_file):
    """
    Read the coverage file and return a dictionary of the coverage information.
    """
    from collections import defaultdict
    cov_dict = defaultdict()
    with open(coverage_file) as f:
        for line in f:
            chrom, pos, cov = line.strip().split()
            cov_dict[str(pos)] = int(cov)
    return cov_dict

def open_file(file_path):
    import gzip
    if file_path.endswith('.gz'):
        file = gzip.open(file_path, 'rt')
    else:
        file = open(file_path, 'r')
    return file

def vcf2dict(vcf_file):
    """
    Read the VCF file and return a dictionary of the VCF information.
    """
    vcf_dict = defaultdict()
    header = '##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n'
    vcff = open_file(vcf_file)
    for line in vcff:
        if line.startswith('#'):
            if not line.startswith('##fileformat'):
                if line.startswith('##FORMAT=<ID=DP,Number'):
                    line = '##FORMAT=<ID=DP,Number=1,Type=Float,Description="Reads covering the REF position">\n'
                    header += line
                else:
                    header += line
        else:
            chrom, pos, id, ref, alt, qual, filter, info, format, sample = line.strip().split('\t')
            vcf_dict[str(pos)] = line.strip()
    vcff.close()
    return vcf_dict, header


def add_wild2vcf(genome, coverage, vcf_file, out_file, count=10):
    """
    fill unkown positions in VCF file with wild type and unkonwn ./. from genome file using coverage information.
    """
    vcf_dict, header = vcf2dict(vcf_file)
    cov_dict = coverage2dict(coverage)
    fa_str = genome2str(genome)
    with open(out_file, 'w') as out:
        out.write(header)
        for i in range(1,len(fa_str)+1):
            if str(i) in vcf_dict:
                print(vcf_dict[str(i)], file=out)
            else:
                chrom = "chrM"
                pos = str(i)
                id = "."
                ref = fa_str[i-1]
                alt = "<NON_REF>"
                qual = "."
                filter = "PASS"
                info = f'END={i}'
                format = "GT:DP:HF:CILOW:CIUP:SDP"
                if cov_dict[str(i)] < count:
                    genome_type = "./."
                else:
                    genome_type = "0/0"
                # if cov_dict[str(i)]:
                DP = cov_dict[str(i)]
                # else:
                #     DP = 0
                sample = f'{genome_type}:{DP}:0:0:0:0'
                print(chrom, pos, id, ref, alt, qual, filter, info, format, sample, sep='\t', file=out)
                
def parse_args():
    example_text= "python vcf2gvcf_addunkonw.py -g mtDNA.fa -c mtDNA.cov -v mtDNA.vcf -o mtDNA_add_wild.vcf -n 10"
    parser = argparse.ArgumentParser(description='Add wild type and unknown ./. to VCF file using coverage information and genome file.',
                                    formatter_class= argparse.RawTextHelpFormatter,
                                    usage = '%(prog)s [-h]',                         
                                    epilog=example_text)
    parser.add_argument('--genome','-g', help='Genome file in FASTA format.')
    parser.add_argument('--coverage','-c', help='Coverage file in BED format.')
    parser.add_argument('--vcf_file','-v',help='VCF file in VCF format.')
    parser.add_argument('--out_file','-o', help='Output VCF file in VCF format.')
    parser.add_argument('--count','-n', type=int, default=10, help='Minimum coverage count to consider a position as known.')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    add_wild2vcf(args.genome, args.coverage, args.vcf_file, args.out_file, args.count)
