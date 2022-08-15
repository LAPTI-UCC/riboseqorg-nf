${task_id} = {"status" : "${status}", "info" : {"script": """${script}""", 
"exit_status": "${exit}",  "workdir": "${workdir}"}}

# nextflow log <run-name> -t scripts/exec_report.py > test.txt
# it's important to keep the model as a .py script and the output as a txt for easiness for interpretation