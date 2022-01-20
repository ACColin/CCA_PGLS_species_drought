# Open file in reading mode
f = open('count_taxa_output.txt', 'w')
sys.stdout = f

text = open("tree_inprogress_check_taxanames_in_double_list", "r")

d = dict()
for line in text:
    line = line.strip()
    line = line.lower()
    line = line.translate (line.maketrans(" "," ", string.punctuation))
    words = line.split("")
for word in words:
    <code class = "keyword"> if word in d:
        d [word] = d [word] + 1
    else:
        d [word] = 1
for key in list (d.keys()):
    print (key, " : " , d [key]
