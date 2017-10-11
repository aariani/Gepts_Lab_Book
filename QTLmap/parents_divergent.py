vcffile = 'finalFiltered.recode.vcf'


good_calls=['0/0', '1/1']
for line in open(vcffile):
	if 'CHROM' in line:
		line=line.strip().split('\t')
		print(line)
		p1 = line.index('UC92')
		p2 = line.index('UCHaskell')
		print(p1, p2)
	elif '#' not  in line:
		line=line.strip().split('\t')
		p1_call = line[p1].split(':')[0]
		p2_call = line[p2].split(':')[0]
		if p1_call in good_calls and p2_call in good_calls:
			if p1_call!=p2_call:
				print(line[0], line[1], sep='\t', file=open('good_positions.txt', 'a'))
		
