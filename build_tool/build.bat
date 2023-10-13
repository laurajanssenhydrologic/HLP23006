:: save working directory
set mypath=%cd% 

:: set correct conda environment in the working directory
call C:\Users\Koen.Reef\Miniconda3\Scripts\activate.bat C:\Users\Koen.Reef\Miniconda3
cd /d %mypath%
call conda activate p1414test

:: convert jupyter-book into sphinx project
call jupyter-book config sphinx ../jupyter_book

:: add paths with python files for autodoc
call python -m add_path

:: clear folders with current build
call RD /S /Q "../docs"
call MD "../docs"
call copy /b NUL "../docs/.nojekyll"

:: build with sphinx
call sphinx-build ../jupyter_book ../docs -b html

:: keep terminal, can be removed
pause(10)

:: launch webpage
start "" ../docs/index.html
