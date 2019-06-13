import sys

for line in (sys.stdin):
    lineSplit = line.split()
    firstFile = lineSplit[1]
    secondFile = lineSplit[3]
    command=""
    if firstFile.endswith("/pandora.consensus.fq.gz"):
        command=f"diff -qs <(zcat {firstFile} | sort) <(zcat {secondFile} | sort)"
    elif firstFile.endswith("/pandora_multisample.matrix") or firstFile.endswith("/pandora_multisample.vcf_ref.fa") or \
         firstFile.endswith("/pandora_multisample_consensus.vcf") or firstFile.endswith("/pandora_multisample_genotyped.vcf"):
        command = f"diff -qs <(sort {firstFile}) <(sort {secondFile})"
    else:
        command = f"diff -qs {firstFile} {secondFile}"

    print(f"echo \"{command}\"")
    print(command)
