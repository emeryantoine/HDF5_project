Voici mon rendu de projet GLCS

les makefile contient les fonctiond necessaires à la compilation des 2 version du projet.

make heat
	pour compiler la premiere version dans plusieurs fichiers
make run
	pour executer cette premiere version
make show
	pour afficher les fichiers hdf5 de sortie

make heat2
	pour compiler la 2eme version avec l'ecriture en parallèle
make run2
	pour executer cette 2eme version
make show2
	pour afficher le contenu du "heat.h5"

Dans le make file, il y a des valeurs width, height et step definit au debut qui peuvent etre changées pour differents usages.


nettoyage :

make clean
	pour supprimer les .o
make clean+
	pour supprimer les .o, executables et les fichiers .h5