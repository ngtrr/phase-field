
# **Phase-field**
<br>

## What is it?

This project is aiming to simulate the mechanism of magnetostriction by LLG equation using phase-field method.

<br>

## How to use?
Before installing the code you need to get **gcc** to compile and **openCV4** to make the png images.<br>

We recommend to use **pkg-config** to link around opencv packages to compile. Please check your link before.

-- 2019/10/7 --<br>
We newly import **fftw3** for calculate dft more fast. please install by terminal before compiling.<br>
To install fftw3 in Ubuntu, type<br>
>`$ sudo apt install libfftw3-3 libfftw3-dev libfftw3-doc`<br>

-- 2019/10/26 --<br>
Importing **Eigen** for calculate the elastic stiffness tensor. please install by terminal before compiling.<br>

<br>
1. Make the directory for clone repository.<br>

>`$ mkdir new_folder_name`<br>
>`$ cd new_folder_name`<br>

2. Clone the repository.<br>

>`$ git clone git@github.com/ngtrr/phase-field.git`<br>

3. Compile and excute it.<br>

>$ g++ -o output magne.cpp `pkg-config --cflags opencv4` `pkg-config --libs opencv4`<br>
>`$ ./output`<br>

## The directory

- [ ] **magne.cpp**   -   main script.<br>
- [x] **martensite.cpp**  -   phase-field simulation about martensite transition using openCV.<br>
- [x] **command.txt**   -   how to compile the code.<br>
- [x] **a.jpg b.jpg c.jpg**   -   images for the initial of simulation.<br>
- [ ] **the_logic.md**   -   describes the logic of this program.<br>