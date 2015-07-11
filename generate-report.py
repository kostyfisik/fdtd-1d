#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2015  Konstantin Ladutenko <kostyfisik@gmail.com>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
import sys,os,shutil,glob,subprocess
path = os.getcwd()
path_fig = os.path.join(path, 'fig')
print(path_fig)
os.chdir(path_fig)
files = []
for fname in glob.iglob('*'):
    files.append(fname)
files.sort()
print(files)

begin_document = """%%\\documentclass[fullscreen=true, handout]{beamer}
\\documentclass[fullscreen=true]{beamer}
\\usetheme[secheader]{Boadilla}
%%\\usetheme{Boadilla}
%%\\usetheme[height=7mm]{Rochester} 
%%\\useinnertheme{rounded}
 \\usecolortheme{seahorse}
%% \\usecolortheme{whale}
%% \\usecolortheme{rose}
%%\\setbeamertemplate{navigation symbols}{\\insertframenavigationsymbol}
%% or
\\beamertemplatenavigationsymbolsempty
%% \\setbeamertemplate[info line]{footline} %%Не работает
%% {\\quad\\strut\\insertsection
%% \\hfill\\insertframenumber/\\inserttotalframenumber\\strut\\quad} 
\\usepackage[utf8]{inputenc}
\\usepackage[T1]{fontenc} %%apt-get install cm-super
%%\\usepackage[english, russian]{babel}
\\usepackage[english]{babel}
\\emergencystretch=35pt %%Чтобы нe было строчек, вылезших на поля
\\usepackage{ragged2e} %% to use \\justifying
\\usepackage{amsmath}
\\usepackage{booktabs}
\\usepackage{multicol}
\\usepackage{listings}
%% \\lstset{language=C++,
%%   backgroundcolor=\\color{Black},
%%    basicstyle=\\color{White}\\tiny\\ttfamily,
%%    keywordstyle=\\color{Orange},
%%    identifierstyle=\\color{Cyan},
%%    stringstyle=\\color{Red}, 
%%    commentstyle=\\color{Green},  
%%    breaklines=true,
%%    breakatwhitespace=true,
%%    tabsize=10,
%%    showstringspaces=false%%
%% }
%%or better
\\lstset{breakatwhitespace,
language=C++,
columns=fullflexible,
keepspaces,
breaklines,
tabsize=2,
backgroundcolor=\\color{white},
basicstyle=\\color{black}\\small\\ttfamily,
keywordstyle=\\color{orange},
identifierstyle=\\color{blue},
stringstyle=\\color{red}, 
commentstyle=\\color{green},  
showstringspaces=false,
extendedchars=true}
\\newcommand{\\transhi}{40} %% Переменная - хорошая прозрачность
\\setbeamercovered{transparent = \\transhi} 
\\setbeamertemplate{caption}[numbered]
\\begin{document}
%% ----------------------------------------------------------------
\\title{Report \\today}
%% ----------------------------------------------------------------
\\begin{frame}{Contents} %%Добавить интриги
\\begin{multicols}{2}
  \\tableofcontents%%[pausesections]
\\end{multicols}
\\end{frame}
\\section{Figures}
"""
def add_figure(fname):
    if fname[-3:] == "txt":
        comment = "%%"
    else:
        comment = ""
    return """
\\subsection{%s}
\\begin{frame}
  \\begin{figure}
    %s\\includegraphics[width=0.75\\textwidth]{fig/%s}%%
    \\caption{%s}
  \\end{figure}
\\end{frame}
"""%(fname, comment,fname,fname)

new_main_file = os.path.join(path, "report-with-figures.tex")
with open(new_main_file, 'w') as file_:
    file_.write(begin_document)            
    for filename in files:
        file_.write(add_figure(filename))
    file_.write("\\end{document}\n")
