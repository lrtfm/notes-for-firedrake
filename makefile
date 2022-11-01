all: 01_firedrake_install.pdf 02_firedrake_intro.pdf

01_firedrake_install.pdf: 01_firedrake_install.ipynb
	jupyter nbconvert --to pdf 01_firedrake_install.ipynb

02_firedrake_intro.pdf: 02_firedrake_intro.ipynb
	jupyter nbconvert --to pdf 02_firedrake_intro.ipynb
