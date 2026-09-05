#!/usr/bin/env bash

set -euo pipefail

readonly VERSION="2.0.0"

print_help() {
    cat <<'EOF'
Calculate theoretical sequencing depth from FASTQ read counts.

Usage:
  read_depth_calcu.sh [options] END_COUNT READ_LENGTH GENOME_SIZE [FASTQ ...]

Arguments:
  END_COUNT    1 for single-end data, 2 for paired-end data
  READ_LENGTH  Read length in bases (a positive integer)
  GENOME_SIZE  Haploid genome size in bases (a positive integer)
  FASTQ        One or more .fastq, .fq, .fastq.gz, or .fq.gz files

Options:
  -o FILE      Output TSV file (default: read_depth.tsv; use - for stdout)
  -h, --help   Show this help message
  -V, --version
               Show the program version

If FASTQ files are omitted, supported FASTQ files in the current directory
are discovered automatically. Each input file produces one output row.

For paired-end data, provide one mate (normally R1) per sample because
END_COUNT=2 already accounts for both ends. Supplying both R1 and R2 will
produce two independent rows for the same library.

Formula:
  depth = read_count * END_COUNT * READ_LENGTH / GENOME_SIZE
EOF
}

die() {
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

require_positive_integer() {
    local name=$1
    local value=$2

    [[ $value =~ ^[1-9][0-9]*$ ]] ||
        die "$name must be a positive integer (received: $value)"
}

count_fastq_records() {
    local file=$1
    local line_count

    if [[ $file == *.gz ]]; then
        if ! line_count=$(gzip -cd -- "$file" | awk 'END { print NR }'); then
            die "could not decompress $file"
        fi
    else
        if ! line_count=$(awk 'END { print NR }' < "$file"); then
            die "could not read $file"
        fi
    fi

    [[ $line_count =~ ^[0-9]+$ ]] || die "could not count lines in $file"
    (( line_count % 4 == 0 )) ||
        die "$file has $line_count lines; a FASTQ file must contain complete four-line records"

    printf '%d\n' "$((line_count / 4))"
}

escape_tsv_field() {
    local value=$1

    value=${value//$'\t'/\\t}
    value=${value//$'\n'/\\n}
    value=${value//$'\r'/\\r}
    printf '%s' "$value"
}

output_file=read_depth.tsv

while (( $# > 0 )); do
    case $1 in
        -o)
            (( $# >= 2 )) || die '-o requires an output file'
            output_file=$2
            shift 2
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        -V|--version)
            printf 'read_depth_calcu.sh %s\n' "$VERSION"
            exit 0
            ;;
        --)
            shift
            break
            ;;
        -*)
            die "unknown option: $1 (use --help for usage)"
            ;;
        *)
            break
            ;;
    esac
done

(( $# >= 3 )) || {
    print_help >&2
    exit 1
}

end_count=$1
read_length=$2
genome_size=$3
shift 3

[[ $end_count == 1 || $end_count == 2 ]] ||
    die "END_COUNT must be 1 or 2 (received: $end_count)"
require_positive_integer READ_LENGTH "$read_length"
require_positive_integer GENOME_SIZE "$genome_size"

files=("$@")
if (( ${#files[@]} == 0 )); then
    shopt -s nullglob
    files=( *.fastq.gz *.fq.gz *.fastq *.fq )
    shopt -u nullglob
fi

(( ${#files[@]} > 0 )) ||
    die 'no FASTQ files found; pass input files explicitly or run in a directory containing FASTQ files'

for file in "${files[@]}"; do
    [[ -f $file ]] || die "input is not a regular file: $file"
    [[ -r $file ]] || die "input is not readable: $file"
    case $file in
        *.fastq|*.fq|*.fastq.gz|*.fq.gz) ;;
        *) die "unsupported input extension: $file" ;;
    esac
done

if [[ $output_file == - ]]; then
    temporary_output=$(mktemp "${TMPDIR:-/tmp}/read-depth.XXXXXX") ||
        die 'could not create a temporary output file'
else
    output_directory=$(dirname -- "$output_file")
    output_name=$(basename -- "$output_file")
    [[ -d $output_directory ]] || die "output directory does not exist: $output_directory"

    output_path=$(cd "$output_directory" && printf '%s/%s' "$PWD" "$output_name")
    for file in "${files[@]}"; do
        input_directory=$(dirname -- "$file")
        input_name=$(basename -- "$file")
        input_path=$(cd "$input_directory" && printf '%s/%s' "$PWD" "$input_name")

        [[ $output_path != "$input_path" ]] ||
            die "output file must not overwrite an input FASTQ: $file"
        if [[ -e $output_file && $output_file -ef $file ]]; then
            die "output file must not refer to an input FASTQ: $file"
        fi
    done

    temporary_output=$(mktemp "$output_directory/.${output_name}.XXXXXX") ||
        die "could not create a temporary file in $output_directory"
fi

cleanup() {
    [[ -n ${temporary_output:-} && -f $temporary_output ]] && rm -f -- "$temporary_output"
    return 0
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

printf 'file\tread_count\tend_count\tread_length_bp\tsequenced_bases\tgenome_size_bp\tdepth_x\n' > "$temporary_output"

for file in "${files[@]}"; do
    read_count=$(count_fastq_records "$file")
    calculations=$(awk \
        -v reads="$read_count" \
        -v ends="$end_count" \
        -v read_len="$read_length" \
        -v genome="$genome_size" \
        'BEGIN {
            bases = reads * ends * read_len
            printf "%.0f\t%.0f\t%.6f", bases, genome, bases / genome
        }')

    escaped_file=$(escape_tsv_field "$file")
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$escaped_file" "$read_count" "$end_count" "$read_length" \
        "$calculations" >> "$temporary_output"
    printf 'Processed %s: %s reads\n' "$file" "$read_count" >&2
done

if [[ $output_file == - ]]; then
    cat "$temporary_output"
else
    mv -f -- "$temporary_output" "$output_file"
    temporary_output=
    printf 'Results written to %s\n' "$output_file" >&2
fi
