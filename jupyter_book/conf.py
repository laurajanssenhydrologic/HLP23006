###############################################################################
# Auto-generated by `jupyter-book config`
# If you wish to continue using _config.yml, make edits to that file and
# re-generate this one.
###############################################################################
author = 'KRG Reef'
comments_config = {'hypothesis': False, 'utterances': False}
copyright = '2022'
exclude_patterns = ['**.ipynb_checkpoints', '.DS_Store', 'Thumbs.db', '_build']
extensions = ['sphinx_togglebutton', 'sphinx_copybutton', 'myst_nb', 'jupyter_book', 'sphinx_thebe', 'sphinx_comments', 'sphinx_external_toc', 'sphinx.ext.intersphinx', 'sphinx_design', 'sphinx_book_theme', 'sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode', 'sphinx.ext.autosummary', 'sphinx_jupyterbook_latex']
external_toc_exclude_missing = False
external_toc_path = '_toc.yml'
html_baseurl = ''
html_favicon = ''
html_logo = 'logo.png'
html_sourcelink_suffix = ''
html_theme = 'sphinx_book_theme'
html_theme_options = {'search_bar_text': 'Search this book...', 'launch_buttons': {'notebook_interface': 'classic', 'binderhub_url': '', 'jupyterhub_url': '', 'thebe': False, 'colab_url': ''}, 'path_to_docs': 'jupyter_book/', 'repository_url': 'https://github.com/HydroLogicBV/P1414', 'repository_branch': 'main', 'extra_footer': '', 'home_page_in_toc': True, 'announcement': '', 'analytics': {'google_analytics_id': ''}, 'use_repository_button': True, 'use_edit_page_button': False, 'use_issues_button': True}
html_title = 'D-HydroLogic'
latex_engine = 'pdflatex'
myst_enable_extensions = ['colon_fence', 'dollarmath', 'linkify', 'substitution', 'tasklist']
myst_url_schemes = ['mailto', 'http', 'https']
nb_execution_allow_errors = False
nb_execution_cache_path = ''
nb_execution_excludepatterns = []
nb_execution_in_temp = False
nb_execution_mode = 'off'
nb_execution_timeout = 30
nb_output_stderr = 'show'
numfig = True
pygments_style = 'sphinx'
suppress_warnings = ['myst.domains']
use_jupyterbook_latex = True
use_multitoc_numbering = True
import os
import sys
sys.path.insert(0, os.path.abspath('..\HydroLogic_Inundation_toolbox\Readers'))
sys.path.insert(0, os.path.abspath('..\HydroLogic_Inundation_toolbox'))
sys.path.insert(0, os.path.abspath('..\Code'))
sys.path.insert(0, os.path.abspath('..\Code\data_structures'))
napoleon_custom_sections = [('Returns', 'params_style')]