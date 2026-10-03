import sys
rows=[]
for line in open(sys.argv[1]):
    p=line.rstrip().split("\t")
    rows.append((p[0],int(p[1])))
lengths=sorted((n for _,n in rows),reverse=True)
total=sum(lengths)
n50=0
acc=0
for n in lengths:
    acc+=n
    if acc>=total/2:
        n50=n
        break
primary={f"Chr{i}" for i in range(1,8)}
pl=sum(n for name,n in rows if name in primary)
print(f"sequences\t{len(rows)}")
print(f"primary_chromosomes_Chr1_to_Chr7\t{sum(name in primary for name,_ in rows)}")
print(f"other_sequences\t{sum(name not in primary for name,_ in rows)}")
print(f"total_bp\t{total}")
print(f"primary_chromosome_bp\t{pl}")
print(f"other_sequence_bp\t{total-pl}")
print(f"N50_bp\t{n50}")
print(f"largest_sequence_bp\t{lengths[0]}")
print("primary_lengths_bp")
for name,n in rows:
    if name in primary:
        print(f"{name}\t{n}")
