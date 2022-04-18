from sys import argv
# python calculating_chrom.sizes.py genome_input.fa output.chrom.sizes


def get_chrom_sizes(genome, output):
	'''
	for the provided genome, create a chromosome sizes file at the given output
	'''
	chromSizesoutput = open(output,"w")

	records = []
	record = False
	for line in open(genome, 'r').readlines():
		if line[0] == '>':
			if record:
				records.append(record)
			record = [line.strip("\n").split(' ')[0][1:], 0]

		else:
			sequence = line.strip('\n')
			record[1] += len(sequence)
			
	for seq_record in records:
		output_line = '%s\t%i\n' % (seq_record[0], seq_record[1])
		chromSizesoutput.write(output_line)

	chromSizesoutput.close()


if __name__ == '__main__':
	genome = str(argv[1])
	output = str(argv[2])

