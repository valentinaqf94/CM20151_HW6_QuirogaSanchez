all: Orbitas-1yr.png Orbitas-1000yr.png orbitas.txt

Orbitas-1yr.png: trayectorias.py

	python trayectorias.py

Orbitas-1000yr.png: trayectorias.py

	python trayectorias.py

trayectorias.py: orbitas.txt

	python trayectorias.py; touch trayctorias.py

orbitas.txt: 4body.x

	./4body.x ic.txt 0.00149597871 1000.0

4body.x: 4body.c

	cc 4body.c -o 4body.x 