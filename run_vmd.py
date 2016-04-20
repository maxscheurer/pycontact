import subprocess
import re


command = 'vmd -dispdev text -e contacts.tcl'
process = subprocess.Popen(command.split(" "), stdout=subprocess.PIPE)
started = False
startedCounting = False
finished = False
i = 1
limit = 50
for line in iter(process.stdout.readline, b''):
    if  not started:
        if re.search(b"starting contacts", line) is not None:
            started = True
    elif started:
        m = re.search(b"[0-9]+", line)
        if m is not None:
            n = m.group(0).decode("utf-8")
            if int(n) == limit and startedCounting and not finished:
                finished = True
                print(i)
                print("Finished")
            if int(n) == i and not finished:
                startedCounting = True
                print(i)
                i += 1
    else:
        print(line.decode("utf-8"))

print("VMD finished")