import os
def combine_text_files(input_files, output_file):
    with open(output_file, 'w') as outfile:
        for fname in input_files:
            equip = fname[0:4]
            with open(os.path.join('output/inv_results/', fname)) as infile:
                for line in infile:
                    if line.strip():
                        outfile.write(line.rstrip('\n') + equip + ',\n')
                outfile.write(infile.read())

files_to_combine = []
# Example usage:
for file in os.listdir('output/inv_results/'):
    if file.endswith('.txt'):
        files_to_combine.append(file)

output_filename = 'combined_python.txt'
combine_text_files(files_to_combine, output_filename)