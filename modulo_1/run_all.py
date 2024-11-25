import subprocess

scripts = [
    "./modulo_1/figure_1.py",
    "./modulo_1/figure_2.py"
    "./modulo_1/figure_3.py",
    "./modulo_1/figure_4_5.py",
    "./modulo_1/figure_6.py",
    "./modulo_1/figure_7.py",
]

for script in scripts:
    print(f"Eseguendo {script}...")
    subprocess.run(["python3", script])
