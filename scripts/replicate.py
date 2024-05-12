import sys
# Description: This python script uses the given data and generates the replicated dataset for future use.
#argv1: file containing the positions and total and annotated genes
#argv2: output file name
#argv3: the name of the individual
#argv4: the group of the individual
#argv5: the tissue from which DNA extracted
#argv6: genetic sex
with open(sys.argv[1], "r") as input:
    with open(sys.argv[2], "w") as output:
        for line in input:
            the_list = line.split()
            count=int(the_list[2])
            if int(the_list[3])>= 4:
                for i in range(int(the_list[3])):
                    output.write(the_list[0] + "\t" + the_list[1] + "\t")
                    if count == 0:
                        output.write("0")
                    else:
                        output.write("1")
                        count-=1
                    output.write("\t" + sys.argv[3] + "\t" + sys.argv[4] + "\t" + the_list[4] + "\t" + sys.argv[5] + "\t" + sys.argv[6]  + "\n")

