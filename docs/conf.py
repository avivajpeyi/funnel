"""Sphinx configuration."""
from datetime import datetime


project = "Funnel"
author = "Avi Vajpeyi"
copyright = f"{datetime.now().year}, {author}"
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
autodoc_typehints = "description"
