from ftplib import FTP
from operator import itemgetter
from sys import argv,stdout

if len(argv) != 2:
    print("ERROR: Incorrect number of arguments")
    print("USAGE: python3 ftp_bacteria.py <bacteria_list_txt>")
    exit(-1)
l = [i.strip() for i in open(argv[1]).read().strip().splitlines()]

print("Connecting to NCBI FTP...", end=' ')
stdout.flush()
ftp = FTP("ftp.ncbi.nlm.nih.gov")
print(ftp.login())
stdout.flush()

print("Changing to RefSeq bacteria directory...", end=' ')
stdout.flush()
ftp.cwd("genomes")
ftp.cwd("refseq")
ftp.cwd("bacteria")
genomes_dir = ftp.pwd()
print("done")
stdout.flush()

failed = set()
print("Getting list of bacteria folders in RefSeq...", end=' ')
stdout.flush()
dirs = [i[0] for i in ftp.mlsd() if i[0] != '.' and i[0] != '..']
print("done")
stdout.flush()
for name in l:
    print("Trying " + name + "...", end=' ')
    stdout.flush()
    query = name.lower().replace(' ','_')
    found = False
    for d in dirs:
        if query in d.lower():
            ftp.cwd(d)
            found = True
            break
    if found:
        ftp.cwd("latest_assembly_versions")
        assemblies = sorted([(i[0], i[1]['modify']) for i in ftp.mlsd() if i[0] != '.' and i[0] != '..'], key=itemgetter(1))
        ftp.cwd(assemblies[0][0])
        filename = assemblies[0][0] + '_genomic.fna.gz'
        f = open(name + '.fna.gz','wb')
        ftp.retrbinary('RETR %s' % filename, f.write)
        f.close()
        print("Success")
        stdout.flush()
        ftp.cwd(genomes_dir)
    else:
        failed.add(name)
        print("Fail")
        stdout.flush()

print()
if len(failed) == 0:
    print("All genomes were downloaded successfully")
else:
    open("failed.txt",'w').write('\n'.join(sorted(failed)))
    print("Genomes that failed to download were written in 'failed.txt'")
