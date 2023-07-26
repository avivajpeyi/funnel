Search.setIndex({"docnames": ["examples/linear_regression", "examples/lvk", "examples/lvk_injection", "index"], "filenames": ["examples/linear_regression.ipynb", "examples/lvk.ipynb", "examples/lvk_injection.ipynb", "index.md"], "titles": ["EX1: Linear regression", "EX3: GW150914", "EX2: CBC-GW Injection", "Funnel"], "terms": {"we": [0, 1, 2], "first": 0, "consid": 0, "simpl": 0, "model": [0, 1, 2], "gaussian": 0, "likelihood": [0, 1, 2, 3], "comput": 0, "evid": [0, 1, 3], "fourier": [0, 3], "integr": [0, 3], "method": [0, 3], "load_ext": [0, 1, 2], "autoreload": [0, 1, 2], "2": [0, 1, 2], "matplotlib": [0, 1, 2], "inlin": [0, 1, 2], "import": [0, 1, 2], "os": [0, 1, 2], "shutil": [0, 1, 2], "numpi": [0, 1, 2], "np": [0, 1, 2], "warn": 0, "from": [0, 1, 2], "pyplot": [0, 2], "plt": [0, 2], "bilbi": [0, 2], "log": [0, 1, 2], "funnel": [0, 1, 2], "plot": [0, 1, 2], "plot_fi_evidence_result": [0, 1, 2], "fi_cor": [0, 1, 2], "get_fi_lnz_list": [0, 1, 2], "clean": [0, 1, 2], "fals": [0, 1, 2], "outdir": [0, 1, 2], "out_lin": 0, "path": [0, 1, 2], "exist": [0, 1, 2], "rmtree": [0, 1, 2], "makedir": [0, 1, 2], "exist_ok": [0, 1, 2], "true": [0, 1, 2], "random": [0, 1, 2], "seed": [0, 1, 2], "42": [0, 1, 2], "filterwarn": 0, "ignor": 0, "getlogg": [0, 2], "setlevel": [0, 2], "critic": [0, 2], "defin": 0, "our": 0, "signal": [0, 2], "thi": [0, 1, 2], "case": 0, "function": 0, "def": [0, 1], "time": 0, "m": 0, "c": 0, "return": [0, 1], "now": 0, "inject": [0, 3], "paramet": [0, 1, 2], "which": 0, "make": 0, "simul": [0, 2], "injection_paramet": [0, 2], "dict": [0, 1, 2], "0": [0, 1, 2], "5": [0, 1, 2], "sampling_frequ": [0, 2], "10": [0, 1, 2], "time_dur": 0, "arang": 0, "1": [0, 1, 2], "n": [0, 3], "len": [0, 1], "sigma": 0, "normal": 0, "01": 0, "quickli": 0, "check": [0, 2], "look": 0, "sensibl": 0, "fig": [0, 1, 2], "ax": 0, "subplot": 0, "o": [0, 3], "label": [0, 2], "r": [0, 1, 2], "set_xlabel": 0, "set_ylabel": [0, 1, 2], "y": 0, "legend": 0, "class": 0, "gaussianlikelihood": 0, "analytical1dlikelihood": 0, "__init__": 0, "self": 0, "x": 0, "func": 0, "none": 0, "kwarg": 0, "super": 0, "log_likelihood": [0, 1, 2], "sum": 0, "residu": 0, "pi": 0, "noise_log_likelihood": [0, 2], "prior": [0, 1, 2], "core": 0, "priordict": 0, "uniform": [0, 2], "line": 0, "res_fn": 0, "f": [0, 1, 2, 3], "_result": [0, 2], "json": [0, 2], "result": [0, 2], "read_in_result": 0, "els": [0, 2], "run_sampl": [0, 2], "sampler": [0, 1, 2], "dynesti": [0, 2], "nlive": 0, "1500": 0, "lnz_err": [0, 1, 2], "log_evid": [0, 2], "log_evidence_err": [0, 2], "print": [0, 1, 2], "2f": [0, 1, 2], "lnbf": [0, 1, 2], "log_bayes_factor": [0, 2], "posterior": [0, 2, 3], "head": [0, 1], "143": 0, "80": 0, "08": 0, "462": 0, "81": [0, 1], "log_prior": [0, 1, 2], "357726": 0, "041932": 0, "460": 0, "450103": 0, "995732": 0, "611756": 0, "560525": 0, "124596": 0, "548062": 0, "373219": 0, "476530": 0, "3": [0, 1, 2], "532804": 0, "427149": 0, "706528": 0, "4": [0, 1, 2], "538333": 0, "322802": 0, "891511": 0, "gaussianlikelihoodnonoiselikelihood": 0, "line_wo_nois": 0, "65": 0, "nan": 0, "622600": 0, "653672": 0, "146": 0, "520300": 0, "441560": 0, "881200": 0, "145": 0, "141573": 0, "518018": 0, "263875": 0, "144": 0, "852101": 0, "505545": 0, "197196": 0, "554968": 0, "516669": 0, "530921": 0, "507963": 0, "plot_corn": [0, 2], "suptitl": 0, "note": [0, 1, 2], "store": [0, 2], "ratio": [0, 1, 2], "column": [0, 1, 2], "when": 0, "compar": 0, "abov": 0, "two": 0, "exampl": [0, 1, 2, 3], "lnz_file": [0, 2], "npz": [0, 1, 2], "r_val": [0, 1, 2], "num_ref_param": [0, 1, 2], "100": [0, 1, 2], "cache_fn": [0, 1, 2], "plt_kwg": [0, 1, 2], "sampling_lnz": [0, 1, 2], "tight_layout": [0, 2], "plot_all_lnz": [0, 1, 2], "gca": [0, 1, 2], "set_xlim": 0, "1e3": 0, "set_ylim": 0, "147": 0, "137": 0, "In": [1, 2], "an": [1, 2], "lvk": [1, 2], "zenodo": [1, 2], "cbc": [1, 3], "download": 1, "h5py": 1, "panda": 1, "pd": 1, "collect": 1, "namedtupl": 1, "fpath": 1, "igwn": 1, "gwtc2p1": 1, "v2": 1, "gw150914_095045_pedatarelease_mixed_cosmo": 1, "h5": 1, "lvk_data": 1, "load_lvk_data": 1, "file": 1, "sampling_param": 1, "list": 1, "c01": 1, "imrphenomxphm": 1, "analyt": 1, "kei": [1, 2], "sampler_data": 1, "meta_data": 1, "ln_evid": 1, "ln_evidence_error": 1, "ln_bayes_factor": 1, "post": [1, 2], "posterior_sampl": 1, "datafram": 1, "name": 1, "dtype": 1, "onli": 1, "keep": 1, "more": 1, "than": [1, 2], "one": 1, "uniqu": 1, "valu": 1, "ie": 1, "drop": 1, "delta": 1, "loc": 1, "nuniqu": 1, "gw150914_data": 1, "shape": 1, "147634": 1, "60": 1, "6984": 1, "67": 1, "14": 1, "303": 1, "45": 1, "a_1": [1, 2], "a_2": [1, 2], "azimuth": 1, "chirp_mass": [1, 2], "geocent_tim": [1, 2], "luminosity_dist": [1, 2], "mass_1": [1, 2], "mass_2": [1, 2], "mass_ratio": [1, 2], "phase": [1, 2], "recalib_l1_phase_7": 1, "recalib_l1_phase_8": 1, "recalib_l1_phase_9": 1, "theta_jn": [1, 2], "tilt_1": [1, 2], "tilt_2": [1, 2], "time_jitt": 1, "zenith": 1, "924210": 1, "331092": 1, "278707": 1, "29": [1, 2], "180724": 1, "126259e": 1, "09": [1, 2], "323": 1, "733223": 1, "37": 1, "815124": 1, "796750": 1, "787959": 1, "840744": 1, "007216": 1, "008044": 1, "012153": 1, "775318": 1, "913304": 1, "676054": 1, "000228": 1, "326161": 1, "322": 1, "304271": 1, "76": 1, "937960": 1, "647369": 1, "313305": 1, "602855": 1, "953047": 1, "510": 1, "955362": 1, "032017": 1, "32": 1, "002137": 1, "864175": 1, "337254": 1, "003571": 1, "001302": 1, "007935": 1, "660558": 1, "841044": 1, "839976": 1, "000325": 1, "334247": 1, "320": 1, "301483": 1, "78": 1, "620833": 1, "205678": 1, "875008": 1, "715020": 1, "31": 1, "433890": 1, "500": 1, "928555": 1, "39": 1, "143076": 1, "33": 1, "351053": 1, "852029": 1, "757881": 1, "002202": 1, "005764": 1, "002742": 1, "493798": 1, "365895": 1, "369656": 1, "000120": 1, "345314": 1, "654009": 1, "031577": 1, "711251": 1, "004800": 1, "380721": 1, "30": 1, "741031": 1, "578": 1, "011396": 1, "35": 1, "664826": 1, "34": 1, "963693": 1, "980341": 1, "723087": 1, "009888": 1, "005924": 1, "010894": 1, "967835": 1, "672429": 1, "627161": 1, "000270": 1, "369513": 1, "324": 1, "265661": 1, "326038": 1, "250640": 1, "227152": 1, "783846": 1, "270597": 1, "587": 1, "119007": 1, "250081": 1, "647430": 1, "930130": 1, "542601": 1, "000942": 1, "014119": 1, "010260": 1, "074509": 1, "373358": 1, "525285": 1, "000435": 1, "339495": 1, "729426": 1, "62": 1, "266597": 1, "row": 1, "actual": 1, "lnl": [1, 2], "nois": [1, 2], "out_gw150914_downsampl": 1, "n_samp": 1, "10000": 1, "n_ref_point": 1, "weight": 1, "exp": 1, "try": 1, "out": 1, "refer": 1, "geomspac": 1, "1e9": 1, "weight_samples_by_lnl": 1, "34mfunnel": [1, 2], "0m": [1, 2], "info": [1, 2], "32mcalcul": [1, 2], "size": [1, 2], "58": 1, "6": [1, 2], "77": 1, "alpha": 1, "out_gw150914": 1, "lw": 1, "load": 2, "The": [2, 3], "extens": 2, "alreadi": 2, "To": 2, "reload": 2, "us": [2, 3], "reload_ext": 2, "out_gwinj": 2, "bilby_logg": 2, "durat": 2, "sampling_freq": 2, "min_freq": 2, "1024": 2, "20": 2, "36": 2, "mass": 2, "phi_12": 2, "phi_jl": 2, "spin": 2, "ra": 2, "375": 2, "dec": 2, "2108": 2, "2000": 2, "7": 2, "extrins": 2, "psi": 2, "659": 2, "1126259642": 2, "413": 2, "inj_m1": 2, "inj_m2": 2, "inj_chirp_mass": 2, "convers": 2, "component_masses_to_chirp_mass": 2, "inj_q": 2, "component_masses_to_mass_ratio": 2, "waveform_gener": 2, "waveformgener": 2, "frequency_domain_source_model": 2, "sourc": 2, "lal_binary_black_hol": 2, "parameter_convers": 2, "convert_to_lal_binary_black_hole_paramet": 2, "waveform_argu": 2, "waveform_approxim": 2, "imrphenomd": 2, "reference_frequ": 2, "minimum_frequ": 2, "detector": 2, "ligo": 2, "hanford": 2, "h1": 2, "design": 2, "sensit": 2, "ifo": 2, "interferometerlist": 2, "set_strain_data_from_power_spectral_dens": 2, "start_tim": 2, "inject_sign": 2, "chirp": 2, "howev": 2, "ar": 2, "quit": 2, "un": 2, "astrophys": 2, "process": 2, "convert": 2, "compon": 2, "bbhpriordict": 2, "uniformincomponentschirpmass": 2, "minimum": 2, "maximum": 2, "perform": 2, "doe": 2, "extend": 2, "space": 2, "longer": 2, "data": 2, "validate_prior": 2, "gravitationalwavetransi": 2, "interferomet": 2, "res_fnam": 2, "cbcresult": 2, "from_json": 2, "filenam": 2, "npoint": 2, "1000": 2, "dlogz": 2, "conversion_funct": 2, "generate_all_bbh_paramet": 2, "result_class": 2, "1984": 2, "70": 2, "updat": 2, "log_likelihood_ratio": 2, "1964": 2, "89": 2, "2054": 2, "98": 2, "90": 2, "figur": 2, "figsiz": 2, "hist": 2, "histtyp": 2, "step": 2, "xlabel": 2, "ytick": 2, "titl": 2, "noiselnl": 2, "2653": 2, "ln": 2, "bf": 2, "u": 3, "rier": 3, "i": 3, "tegratio": 3, "e": 3, "videnc": 3, "calcu": 3, "l": 3, "ation": 3, "python": 3, "implement": 3, "fi": 3, "marginalis": 3, "approxim": 3, "describ": 3, "rotiroti": 3, "et": 3, "al": 3, "2018": 3, "follow": 3, "demonstr": 3, "how": 3, "calcul": 3, "sampl": 3, "three": 3, "differ": 3, "scenario": 3, "linear": 3, "regress": 3, "gw": 3, "gw150914": 3}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"ex1": 0, "linear": 0, "regress": 0, "gener": 0, "some": 0, "data": 0, "nest": [0, 1, 2], "sampl": [0, 1, 2], "lnz": [0, 1, 2], "calcul": [0, 2], "fi": [0, 1, 2], "ex3": 1, "gw150914": 1, "load": 1, "posterior": 1, "comput": [1, 2], "us": 1, "downsampl": 1, "speed": 1, "full": 1, "few": 1, "refernc": 1, "point": 1, "ex2": 2, "cbc": 2, "gw": 2, "inject": 2, "evid": 2, "funnel": 3}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})