# samclipy
Filter SAM file for soft-clipped alignments

## Examples

```bash
# Read in.sam and write output to stout
python samclipy.py in.sam

# Read input from stdin and write stdout to file
python samclipy.py < in.sam > out.sam

# Filter alignments from BAM file, write sorted output to file
samtools view -h in.bam | samclipy.py | samtools sort > out.bam

# Filter alignments from BAM file, write all soft-clipped alignments to file
samtools view -h in.bam | samclipy.py --invert | samtools sort > out.bam
```

## Dependencies

- `re`
- `argparse`

## Full usage

```
python samclipy.py -h
usage: samclipy [-h] [--minClip MINCLIP] [--invert] [samfile]

Filter SAM file for clipped alignments .

positional arguments:
  samfile

optional arguments:
  -h, --help         show this help message and exit
  --minClip MINCLIP  Required total (left + right) clip length.
  --invert           Output only soft-clipped SAM records and ignore the good
                     ones.
```

## Heavily inspired by

- [teloclip](https://github.com/Adamtaranto/teloclip)
  - `splitCIGAR`, `checkClips`, `lenCIGAR`
- [samclip](https://github.com/tseemann/samclip)