#!/usr/bin/env python
# -*- coding: UTF-8 -*-


## install requests  pip install requests -i https://mirrors.aliyun.com/pypi/simple/

genome_ver = "hg38"
region     = "chrM:4321-5678"

c_name = region.split(":")[0]
k1 = region.split(":")[1]

start = k1.split('-')[0]
end   = k1.split('-')[1]

import requests

api_url = "https://api.genome.ucsc.edu/getData/sequence?genome=" + genome_ver +";chrom=" + c_name + ";start=" + start + ";end=" + end

## region     = "genome=hg38;chrom=chrM;start=4321;end=5678"

response = requests.get( api_url )
DNA_seq = response.json()["dna"].upper()

print(DNA_seq)
