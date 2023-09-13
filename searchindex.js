Search.setIndex({"docnames": ["examples/GW150914", "examples/lartillot", "examples/linear_regression", "examples/lnz_comparison", "examples/lvk_injection", "index"], "filenames": ["examples/GW150914.ipynb", "examples/lartillot.ipynb", "examples/linear_regression.ipynb", "examples/lnz_comparison.ipynb", "examples/lvk_injection.ipynb", "index.md"], "titles": ["EX3: GW150914", "EX2: Lartillot\u2019s Model", "EX1: Linear regression", "LnZ comparison", "EX2: CBC-GW Injection", "Funnel"], "terms": {"In": [0, 4], "thi": [0, 1, 4], "exampl": [0, 4, 5], "we": [0, 1, 2, 4], "an": [0, 3, 4], "lvk": [0, 4], "from": [0, 2, 3, 4], "zenodo": [0, 4], "evid": [0, 2, 3, 5], "cbc": [0, 5], "model": [0, 2, 4], "download": [0, 4], "load_ext": [0, 2, 4], "autoreload": [0, 2, 4], "2": [0, 1, 2, 3], "matplotlib": [0, 2, 3, 4], "inlin": [0, 2, 4], "import": [0, 2, 3, 4], "h5py": 0, "numpi": [0, 2, 3, 4], "np": [0, 2, 3, 4], "panda": [0, 3], "pd": [0, 3], "collect": [0, 4], "namedtupl": [0, 4], "funnel": [0, 2, 3, 4], "plot": [0, 1, 2, 4], "plot_fi_evidence_result": [0, 2, 4], "fi_cor": [0, 2, 3, 4], "get_fi_lnz_list": [0, 2, 3, 4], "shutil": [0, 2, 4], "os": [0, 2, 4], "random": [0, 2, 3, 4], "seed": [0, 1, 2, 3, 4], "42": [0, 1, 2, 3, 4], "fpath": [0, 4], "igwn": 0, "gwtc2p1": 0, "v2": 0, "gw150914_095045_pedatarelease_mixed_cosmo": 0, "h5": 0, "lvk_data": [0, 4], "lnz_err": [0, 2, 4], "lnbf": [0, 2, 4], "def": [0, 2, 3, 4], "load_lvk_data": [0, 4], "file": [0, 3], "r": [0, 1, 2, 4], "f": [0, 2, 3, 4, 5], "sampling_param": [0, 4], "list": [0, 1, 3], "c01": 0, "imrphenomxphm": 0, "prior": [0, 1, 2, 3, 4], "analyt": 0, "kei": [0, 3, 4], "remov": [0, 2, 3], "mass_1": [0, 4], "mass_2": [0, 4], "alreadi": [0, 2], "have": 0, "chirp_mass": [0, 4], "mass_ratio": [0, 4], "p": [0, 1], "sampler_data": 0, "meta_data": 0, "sampler": [0, 2, 3, 4], "ln_evid": [0, 3], "0": [0, 1, 2, 3, 4], "ln_evidence_error": 0, "ln_bayes_factor": [0, 3], "post": [0, 1, 3, 4], "posterior_sampl": 0, "datafram": [0, 3], "name": [0, 3], "dtype": 0, "log_likelihood": [0, 2, 3, 4], "log_prior": [0, 2, 3, 4], "onli": [0, 3], "keep": [0, 3], "paramet": [0, 1, 3], "more": 0, "than": [0, 3, 4], "one": 0, "uniqu": 0, "valu": [0, 1, 3], "ie": 0, "drop": 0, "delta": 0, "loc": [0, 2, 4], "nuniqu": [0, 4], "1": [0, 1, 2, 3, 4], "return": [0, 1, 2, 3, 4], "gw150914_data": 0, "get_config": 0, "configs_dict": 0, "config": 0, "k": [0, 3], "v": [0, 1, 3], "item": [0, 3], "str": [0, 3], "margin": [0, 1], "sort": 0, "df": 0, "index": 0, "val": 0, "t": [0, 1, 2, 3, 4], "distance_margin": 0, "b": 0, "true": [0, 1, 2, 3, 4], "phase_margin": 0, "fals": [0, 2, 3, 4], "time_margin": 0, "gw": [0, 5], "150914_data": 0, "column": [0, 2, 4], "without": 0, "set": [0, 1, 2, 3], "cbc_param": 0, "print": [0, 1, 2, 3, 4], "param": [0, 3], "count": 0, "i": [0, 1, 2, 3, 5], "calib": 0, "append": [0, 3], "calib_param": 0, "a_1": [0, 4], "a_2": [0, 4], "3": [0, 1, 2, 3, 4], "azimuth": [0, 4], "4": [0, 1, 2, 3, 4], "5": [0, 2, 3, 4], "geocent_tim": [0, 4], "6": [0, 2, 3, 4], "luminosity_dist": [0, 4], "7": [0, 1, 2, 3, 4], "8": [0, 3], "phase": [0, 4], "9": [0, 2, 3], "phi_12": [0, 4], "10": [0, 2, 3, 4], "phi_jl": [0, 4], "11": [0, 3, 4], "psi": [0, 4], "12": [0, 2, 3], "theta_jn": [0, 4], "13": [0, 3], "tilt_1": [0, 4], "14": [0, 3], "tilt_2": [0, 4], "15": [0, 3], "time_jitt": [0, 4], "16": [0, 3], "zenith": [0, 4], "recalib_h1_amplitude_0": 0, "recalib_h1_amplitude_1": 0, "recalib_h1_amplitude_2": 0, "recalib_h1_amplitude_3": 0, "recalib_h1_amplitude_4": 0, "recalib_h1_amplitude_5": 0, "recalib_h1_amplitude_6": 0, "recalib_h1_amplitude_7": 0, "recalib_h1_amplitude_8": 0, "recalib_h1_amplitude_9": 0, "recalib_h1_phase_0": 0, "recalib_h1_phase_1": 0, "recalib_h1_phase_2": 0, "recalib_h1_phase_3": 0, "recalib_h1_phase_4": 0, "recalib_h1_phase_5": 0, "17": [0, 3, 4], "recalib_h1_phase_6": 0, "18": [0, 3], "recalib_h1_phase_7": 0, "19": [0, 3], "recalib_h1_phase_8": 0, "20": [0, 3, 4], "recalib_h1_phase_9": 0, "21": [0, 3], "recalib_l1_amplitude_0": 0, "22": [0, 3], "recalib_l1_amplitude_1": 0, "23": [0, 3], "recalib_l1_amplitude_2": 0, "24": [0, 3], "recalib_l1_amplitude_3": 0, "25": [0, 1, 3], "recalib_l1_amplitude_4": 0, "26": [0, 2, 3, 4], "recalib_l1_amplitude_5": 0, "27": [0, 3], "recalib_l1_amplitude_6": 0, "28": [0, 3], "recalib_l1_amplitude_7": 0, "29": [0, 3, 4], "recalib_l1_amplitude_8": 0, "30": [0, 3], "recalib_l1_amplitude_9": 0, "31": [0, 3], "recalib_l1_phase_0": 0, "32": [0, 3], "recalib_l1_phase_1": 0, "33": [0, 3], "recalib_l1_phase_2": 0, "34": [0, 3], "recalib_l1_phase_3": 0, "35": [0, 3], "recalib_l1_phase_4": 0, "36": [0, 3, 4], "recalib_l1_phase_5": 0, "37": [0, 3], "recalib_l1_phase_6": 0, "38": [0, 3], "recalib_l1_phase_7": 0, "39": [0, 3], "recalib_l1_phase_8": 0, "40": [0, 1, 3], "recalib_l1_phase_9": 0, "data": [0, 1, 4], "skip": 0, "deltafunct": 0, "2d": 0, "els": [0, 2, 3, 4], "close": 0, "uniform": [0, 2, 4], "minimum": [0, 3, 4], "maximum": [0, 3, 4], "99": 0, "latex_label": [0, 3], "unit": [0, 3], "none": [0, 2, 3, 4], "boundari": [0, 3], "283185307179586": 0, "epsilon": 0, "period": [0, 3], "uniformincomponentschirpmass": [0, 4], "418182160215295": 0, "41": [0, 3], "97447913941358": 0, "mathcal": [0, 2], "m": [0, 2], "m_": 0, "odot": 0, "1126259462": 0, "2910001": 0, "491": 0, "t_c": 0, "s": [0, 2, 3], "powerlaw": 0, "alpha": [0, 2, 4], "10000": [0, 2], "d_l": 0, "mpc": 0, "constraint": 0, "1000": [0, 3, 4], "m_1": 0, "m_2": 0, "uniformincomponentsmassratio": 0, "05": [0, 3], "q": 0, "phi": 0, "phi_": 0, "jl": 0, "141592653589793": 0, "gaussian": [0, 2], "mu": [0, 3], "0022497139999999582": 0, "sigma": [0, 1, 2, 3], "032482292999999995": 0, "A": [0, 3], "h1_0": 0, "reflect": [0, 3], "0022525476996070083": 0, "03909717651863342": 0, "h1_1": 0, "004632777733005889": 0, "04078174073734512": 0, "h1_2": 0, "003987700439142069": 0, "03287682593561539": 0, "h1_3": 0, "0018107999614999065": 0, "01936965624646965": 0, "h1_4": 0, "003528123018086278": 0, "012591404744135145": 0, "h1_5": 0, "0018296964612414133": 0, "01075481852283716": 0, "h1_6": 0, "005059244326670863": 0, "008407555954093059": 0, "h1_7": 0, "0010763014001435082": 0, "01391150436119593": 0, "h1_8": 0, "0017064540000000237": 0, "016027046000000045": 0, "h1_9": 0, "009355752999999998": 0, "030123145": 0, "0050309026057799855": 0, "04177354702981138": 0, "0043774230879091115": 0, "045027737962814354": 0, "0023545243645386326": 0, "036339899119904644": 0, "00444209518471385": 0, "032837260835226936": 0, "006872683631006688": 0, "02806419946997496": 0, "01508554971905464": 0, "022470973360749408": 0, "009105656795490394": 0, "0186470280190447": 0, "01327344249855343": 0, "019250972190353967": 0, "006712540999999999": 0, "020072625": 0, "005476730999999901": 0, "08406088199999993": 0, "l1_0": 0, "012104387120680028": 0, "07651583147867605": 0, "l1_1": 0, "013219120050478253": 0, "06508029899775643": 0, "l1_2": 0, "006657329180423924": 0, "05665180922845485": 0, "l1_3": 0, "0026286630232402308": 0, "05041154190241824": 0, "l1_4": 0, "008989482013408759": 0, "0391433695986004": 0, "l1_5": 0, "0112041341409286": 0, "02176964617289734": 0, "l1_6": 0, "007323119754866439": 0, "01662038622968457": 0, "l1_7": 0, "007118273569718694": 0, "017534092484121667": 0, "l1_8": 0, "43": [0, 3], "00654931000000003": 0, "018560549000000037": 0, "l1_9": 0, "44": [0, 3], "027121368999999992": 0, "065822439": 0, "45": [0, 3, 4], "015666055283929895": 0, "06971337917415898": 0, "46": [0, 3], "0021521761445731018": 0, "06846083085274324": 0, "47": [0, 3], "0022536976421022887": 0, "06194252192559962": 0, "48": [0, 3], "002461783883218465": 0, "048859883075223744": 0, "49": [0, 3], "0067988571593666286": 0, "03330895148538794": 0, "50": [0, 3], "0015671495057928866": 0, "0194164413063868": 0, "51": [0, 3], "0016729416638600334": 0, "013796249113051186": 0, "52": [0, 3], "003950310144720625": 0, "014338015652101466": 0, "53": [0, 3], "008061858000000002": 0, "015019266": 0, "54": [0, 3], "sine": 0, "theta_": 0, "jn": 0, "55": [0, 3], "theta_1": [0, 2], "56": [0, 3, 4], "theta_2": [0, 2], "57": [0, 3], "00048828125": 0, "58": [0, 1, 3], "kappa": 0, "recalib_h1_frequency_0": 0, "peak": 0, "999999999999996": 0, "recalib_h1_frequency_1": 0, "514434824844255": 0, "recalib_h1_frequency_2": 0, "55653663398341": 0, "recalib_h1_frequency_3": 0, "71": 0, "03232013940804": 0, "recalib_h1_frequency_4": 0, "108": 0, "37555516757195": 0, "recalib_h1_frequency_5": 0, "165": 0, "35094073835938": 0, "recalib_h1_frequency_6": 0, "252": 0, "27952521936786": 0, "recalib_h1_frequency_7": 0, "384": 0, "9083564974527": 0, "recalib_h1_frequency_8": 0, "587": 0, "2630478939719": 0, "recalib_h1_frequency_9": 0, "895": 0, "9999999999999": 0, "recalib_l1_frequency_0": 0, "recalib_l1_frequency_1": 0, "recalib_l1_frequency_2": 0, "recalib_l1_frequency_3": 0, "recalib_l1_frequency_4": 0, "recalib_l1_frequency_5": 0, "recalib_l1_frequency_6": 0, "recalib_l1_frequency_7": 0, "recalib_l1_frequency_8": 0, "recalib_l1_frequency_9": 0, "shape": [0, 4], "2f": [0, 2, 3, 4], "head": [0, 2, 4], "147634": 0, "6984": 0, "67": 0, "303": 0, "242101e": 0, "01": [0, 1, 2, 3, 4], "473688e": 0, "056780e": 0, "112512e": 0, "506404e": 0, "310924e": 0, "133050e": 0, "750084e": 0, "800071e": 0, "03": [0, 3], "271518e": 0, "278707e": 0, "00": [0, 3, 4], "602855e": 0, "715020e": 0, "380721e": 0, "783846e": 0, "918072e": 0, "995305e": 0, "143389e": 0, "074103e": 0, "127060e": 0, "126259e": 0, "09": [0, 3], "237332e": 0, "02": [0, 3], "109554e": 0, "009286e": 0, "780114e": 0, "871190e": 0, "879585e": 0, "641748e": 0, "520294e": 0, "803410e": 0, "301303e": 0, "840744e": 0, "337254e": 0, "757881e": 0, "723087e": 0, "542601e": 0, "394340e": 0, "868170e": 0, "916934e": 0, "591901e": 0, "435515e": 0, "210051e": 0, "486261e": 0, "576937e": 0, "338661e": 0, "019726e": 0, "800080e": 0, "571376e": 0, "272124e": 0, "861790e": 0, "418739e": 0, "999532e": 0, "472740e": 0, "900000e": 0, "081133e": 0, "977923e": 0, "383615e": 0, "474793e": 0, "200648e": 0, "799967e": 0, "687969e": 0, "688680e": 0, "209578e": 0, "758393e": 0, "357331e": 0, "604320e": 0, "698306e": 0, "626471e": 0, "131780e": 0, "152776e": 0, "927244e": 0, "305150e": 0, "121430e": 0, "399314e": 0, "998931e": 0, "089422e": 0, "041142e": 0, "582923e": 0, "287890e": 0, "776959e": 0, "428993e": 0, "394128e": 0, "04": [0, 3], "274103e": 0, "444985e": 0, "033461e": 0, "928833e": 0, "628103e": 0, "094559e": 0, "589106e": 0, "269532e": 0, "055759e": 0, "309187e": 0, "280925e": 0, "110290e": 0, "604852e": 0, "810975e": 0, "221316e": 0, "607438e": 0, "675564e": 0, "992964e": 0, "210042e": 0, "719320e": 0, "060107e": 0, "166049e": 0, "323885e": 0, "256244e": 0, "345580e": 0, "485112e": 0, "924302e": 0, "663208e": 0, "985714e": 0, "148825e": 0, "037495e": 0, "528264e": 0, "690566e": 0, "084741e": 0, "833606e": 0, "235240e": 0, "403403e": 0, "949062e": 0, "622115e": 0, "974057e": 0, "740506e": 0, "063412e": 0, "390864e": 0, "773227e": 0, "427238e": 0, "936114e": 0, "498003e": 0, "420854e": 0, "673151e": 0, "205619e": 0, "312319e": 0, "284825e": 0, "316393e": 0, "395785e": 0, "259316e": 0, "696564e": 0, "007685e": 0, "100947e": 0, "684205e": 0, "409629e": 0, "934225e": 0, "075666e": 0, "741046e": 0, "409607e": 0, "438222e": 0, "741232e": 0, "067432e": 0, "260970e": 0, "605238e": 0, "149459e": 0, "533687e": 0, "256302e": 0, "926279e": 0, "115813e": 0, "511059e": 0, "943587e": 0, "095780e": 0, "852026e": 0, "819107e": 0, "530716e": 0, "390149e": 0, "031972e": 0, "972729e": 0, "843934e": 0, "488562e": 0, "857405e": 0, "939737e": 0, "916733e": 0, "209214e": 0, "366417e": 0, "900651e": 0, "958842e": 0, "124533e": 0, "328298e": 0, "566428e": 0, "372492e": 0, "438995e": 0, "581318e": 0, "917951e": 0, "862231e": 0, "591929e": 0, "370625e": 0, "947175e": 0, "395698e": 0, "773761e": 0, "311472e": 0, "900587e": 0, "826212e": 0, "684315e": 0, "266035e": 0, "677954e": 0, "511803e": 0, "485830e": 0, "241593e": 0, "614393e": 0, "712452e": 0, "511530e": 0, "721168e": 0, "419746e": 0, "737868e": 0, "073929e": 0, "091229e": 0, "545232e": 0, "465972e": 0, "689561e": 0, "122678e": 0, "452316e": 0, "197389e": 0, "613542e": 0, "519969e": 0, "114218e": 0, "653476e": 0, "637367e": 0, "982944e": 0, "189141e": 0, "536081e": 0, "243102e": 0, "364353e": 0, "765010e": 0, "976252e": 0, "251552e": 0, "996706e": 0, "091622e": 0, "611963e": 0, "161828e": 0, "238076e": 0, "706364e": 0, "364514e": 0, "986921e": 0, "086198e": 0, "332725e": 0, "956626e": 0, "210257e": 0, "876330e": 0, "216118e": 0, "571173e": 0, "202440e": 0, "887829e": 0, "418501e": 0, "043966e": 0, "301929e": 0, "764206e": 0, "924480e": 0, "411924e": 0, "215333e": 0, "934553e": 0, "742098e": 0, "089411e": 0, "026025e": 0, "775318e": 0, "660558e": 0, "493798e": 0, "967835e": 0, "074509e": 0, "913304e": 0, "841044e": 0, "365895e": 0, "672429e": 0, "373358e": 0, "676054e": 0, "839976e": 0, "369656e": 0, "271606e": 0, "525285e": 0, "283661e": 0, "247132e": 0, "200043e": 0, "701632e": 0, "351720e": 0, "326161e": 0, "334247e": 0, "345314e": 0, "369513e": 0, "339495e": 0, "223043e": 0, "203015e": 0, "206540e": 0, "242657e": 0, "227294e": 0, "693796e": 0, "862083e": 0, "603158e": 0, "132604e": 0, "226660e": 0, "note": [0, 4], "log": [0, 1, 2, 3, 4], "likelihood": [0, 1, 2, 3, 4, 5], "actual": [0, 4], "lnl": [0, 3, 4], "nois": [0, 3, 4], "ratio": [0, 4], "run_fi_comput": [0, 4], "cbc_data": [0, 4], "frac_samp": [0, 4], "n_ref_point": [0, 4], "outdir": [0, 2, 3, 4], "out_inj_downsampl": [0, 4], "clean": [0, 2, 3, 4], "path": [0, 2, 4], "exist": [0, 2, 4], "rmtree": [0, 2, 4], "makedir": [0, 2, 4], "exist_ok": [0, 2, 4], "n_total": [0, 4], "len": [0, 2, 4], "n_samp": [0, 4], "int": [0, 3, 4], "weight": [0, 4], "exp": [0, 1, 4], "100": [0, 2, 3, 4], "refer": [0, 2], "r_val": [0, 2, 3, 4], "_": [0, 2, 4], "num_ref_param": [0, 2, 3, 4], "geomspac": [0, 4], "1e4": [0, 4], "cache_fn": [0, 2, 4], "npz": [0, 2, 4], "weight_samples_by_lnl": [0, 2, 3, 4], "parm": 0, "copi": 0, "out_gw150914_downsampl": 0, "plt_kwg": [0, 2, 4], "dict": [0, 2, 3, 4], "sampling_lnz": [0, 2, 4], "fig": [0, 2, 4], "plot_all_lnz": [0, 2, 4], "gca": [0, 2, 4], "set_ylabel": [0, 2, 4], "34mfunnel": [0, 2, 4], "0m": [0, 2, 4], "info": [0, 2, 3, 4], "32mcalcul": [0, 2, 4], "size": [0, 2, 4], "14763": 0, "32mposterior": [0, 2], "old": 0, "out_gw150914_downsampled_no_calib": 0, "max": [0, 2], "try": 0, "out": 0, "fi_samp": [0, 2], "1e2": 0, "1e5": 0, "77": 0, "plot_corner_and_mark_sampl": 0, "samp": [0, 3], "out_gw150914": 0, "keyboardinterrupt": 0, "traceback": 0, "most": 0, "recent": 0, "call": 0, "last": [0, 3], "cell": 0, "line": [0, 1, 4], "document": 0, "project": 0, "src": 0, "py": [0, 2, 3], "142": 0, "135": [0, 2], "post_mask": 0, "get_post_mask": 0, "refi": 0, "136": 0, "fi_kwarg": 0, "137": 0, "138": 0, "ref_samp": 0, "139": 0, "ref_lnpri": 0, "ln_pri": 0, "140": [0, 2], "ref_lnl": 0, "ln_lnl": 0, "141": 0, "arrai": [0, 3], "fi_ln_evid": 0, "ri": 0, "143": [0, 2], "median_lnz": 0, "nanmedian": 0, "146": 0, "pbar": 0, "set_postfix_str": 0, "med_": 0, "listcomp": 0, "62": [0, 2], "ndarrai": 0, "float": 0, "approx": 0, "some": 0, "The": [0, 1, 2, 5], "approxim": [0, 5], "base": [0, 3], "densiti": [0, 1, 4], "estim": [0, 1, 3], "method": [0, 2, 3, 5], "describ": [0, 5], "60": [0, 3], "61": 0, "approx_ln_post": 0, "fi_ln_posterior": 0, "63": [0, 4], "reference_sampl": 0, "patricio": 0, "suggest": 0, "differ": [0, 1, 2, 5], "each": 0, "max_diff": 0, "nanmax": 0, "diff_from_ref": 0, "axi": 0, "pi": [0, 1, 2], "now": 0, "vector": 0, "sin_diff": 0, "sin": [0, 1], "sum_r": 0, "nansum": 0, "nanprod": 0, "n_dim": [0, 3], "const": 0, "power": [0, 3], "lw": [0, 4], "equat": 1, "theta": 1, "d": [1, 3], "frac": 1, "instal": 1, "packag": [1, 3], "matrixstat": 1, "mvtnorm": 1, "home": [1, 3], "avaj040": [1, 3], "x86_64": 1, "pc": 1, "linux": 1, "gnu": 1, "librari": 1, "lib": [1, 3], "unspecifi": 1, "inflat": 1, "n": [1, 2, 3, 5], "1e3": [1, 2, 4], "rt": 1, "2000": [1, 3, 4], "re": [1, 3], "360": 1, "tau": 1, "eta": 1, "75": 1, "target": 1, "mean": 1, "rep": 1, "var": 1, "diag": 1, "x": [1, 2, 3], "ncol": 1, "nrow": 1, "c": [1, 2], "function": [1, 4], "formula": 1, "given": 1, "right": 1, "after": 1, "ltrue": 1, "delet": 1, "later": 1, "xi": 1, "epanechnikov": 1, "result": [1, 2, 3, 4], "triangl": 1, "doubleexp": 1, "norm": 1, "simul": [1, 4], "na": 1, "300": 1, "ptm": 1, "proc": 1, "time": [1, 2, 3, 4], "joint": 1, "g": 1, "like": 1, "sum": [1, 2, 3], "dnorm": 1, "choos": 1, "evalu": [1, 3], "posterior": [1, 2, 3, 5], "fi": [1, 3, 5], "sampl": [1, 3, 5], "rnorm": 1, "sd": 1, "sqrt": [1, 3], "den": 1, "lpriorlik": 1, "normal": [1, 2, 3], "doubl": 1, "exponenti": 1, "dcauchi": 1, "locat": 1, "scale": [1, 3], "kernel": 1, "co": 1, "matrix": [1, 3], "byrow": 1, "ab": [1, 4], "rowprod": 1, "multipli": 1, "togeth": 1, "dimens": 1, "dmvnorm": 1, "paste0": 1, "iter": 1, "squar": 1, "diff": 1, "png": 1, "gficase4_r": 1, "xlab": 1, "main": 1, "fourier": [1, 2, 5], "integr": [1, 2, 5], "ablin": 1, "col": 1, "red": 1, "mtext": 1, "substitut": 1, "past": 1, "side": 1, "blue": [1, 4], "legend": [1, 2], "topright": 1, "black": 1, "lwd": 1, "lty": 1, "cauchi": 1, "triangular": 1, "dev": 1, "off": 1, "iteration1": 1, "iteration2": 1, "iteration3": 1, "iteration4": 1, "iteration5": 1, "iteration6": 1, "iteration7": 1, "iteration8": 1, "iteration9": 1, "iteration10": 1, "iteration11": 1, "iteration12": 1, "iteration13": 1, "iteration14": 1, "iteration15": 1, "iteration16": 1, "iteration17": 1, "iteration18": 1, "iteration19": 1, "iteration20": 1, "iteration21": 1, "iteration22": 1, "iteration23": 1, "iteration24": 1, "iteration25": 1, "iteration26": 1, "iteration27": 1, "iteration28": 1, "iteration29": 1, "iteration30": 1, "iteration31": 1, "iteration32": 1, "iteration33": 1, "iteration34": 1, "iteration35": 1, "iteration36": 1, "iteration37": 1, "iteration38": 1, "iteration39": 1, "iteration40": 1, "iteration41": 1, "iteration42": 1, "iteration43": 1, "iteration44": 1, "iteration45": 1, "iteration46": 1, "iteration47": 1, "iteration48": 1, "iteration49": 1, "iteration50": 1, "iteration51": 1, "iteration52": 1, "iteration53": 1, "iteration54": 1, "iteration55": 1, "iteration56": 1, "iteration57": 1, "iteration58": 1, "iteration59": 1, "iteration60": 1, "iteration61": 1, "iteration62": 1, "iteration63": 1, "iteration64": 1, "iteration65": 1, "iteration66": 1, "iteration67": 1, "iteration68": 1, "iteration69": 1, "iteration70": 1, "iteration71": 1, "iteration72": 1, "iteration73": 1, "iteration74": 1, "iteration75": 1, "iteration76": 1, "iteration77": 1, "iteration78": 1, "iteration79": 1, "iteration80": 1, "iteration81": 1, "iteration82": 1, "iteration83": 1, "iteration84": 1, "iteration85": 1, "iteration86": 1, "iteration87": 1, "iteration88": 1, "iteration89": 1, "iteration90": 1, "iteration91": 1, "iteration92": 1, "iteration93": 1, "iteration94": 1, "iteration95": 1, "iteration96": 1, "iteration97": 1, "iteration98": 1, "iteration99": 1, "iteration100": 1, "iteration101": 1, "iteration102": 1, "iteration103": 1, "iteration104": 1, "iteration105": 1, "iteration106": 1, "iteration107": 1, "iteration108": 1, "iteration109": 1, "iteration110": 1, "iteration111": 1, "iteration112": 1, "iteration113": 1, "iteration114": 1, "iteration115": 1, "iteration116": 1, "iteration117": 1, "iteration118": 1, "iteration119": 1, "iteration120": 1, "iteration121": 1, "iteration122": 1, "iteration123": 1, "iteration124": 1, "iteration125": 1, "iteration126": 1, "iteration127": 1, "iteration128": 1, "iteration129": 1, "iteration130": 1, "iteration131": 1, "iteration132": 1, "iteration133": 1, "iteration134": 1, "iteration135": 1, "iteration136": 1, "iteration137": 1, "iteration138": 1, "iteration139": 1, "iteration140": 1, "iteration141": 1, "iteration142": 1, "iteration143": 1, "iteration144": 1, "iteration145": 1, "iteration146": 1, "iteration147": 1, "iteration148": 1, "iteration149": 1, "iteration150": 1, "iteration151": 1, "iteration152": 1, "iteration153": 1, "iteration154": 1, "iteration155": 1, "iteration156": 1, "iteration157": 1, "iteration158": 1, "iteration159": 1, "iteration160": 1, "iteration161": 1, "iteration162": 1, "iteration163": 1, "iteration164": 1, "iteration165": 1, "iteration166": 1, "iteration167": 1, "iteration168": 1, "iteration169": 1, "iteration170": 1, "iteration171": 1, "iteration172": 1, "iteration173": 1, "iteration174": 1, "iteration175": 1, "iteration176": 1, "iteration177": 1, "iteration178": 1, "iteration179": 1, "iteration180": 1, "iteration181": 1, "iteration182": 1, "iteration183": 1, "iteration184": 1, "iteration185": 1, "iteration186": 1, "iteration187": 1, "iteration188": 1, "iteration189": 1, "iteration190": 1, "iteration191": 1, "iteration192": 1, "iteration193": 1, "iteration194": 1, "iteration195": 1, "iteration196": 1, "iteration197": 1, "iteration198": 1, "iteration199": 1, "iteration200": 1, "iteration201": 1, "iteration202": 1, "iteration203": 1, "iteration204": 1, "iteration205": 1, "iteration206": 1, "iteration207": 1, "iteration208": 1, "iteration209": 1, "iteration210": 1, "iteration211": 1, "iteration212": 1, "iteration213": 1, "iteration214": 1, "iteration215": 1, "iteration216": 1, "iteration217": 1, "iteration218": 1, "iteration219": 1, "iteration220": 1, "iteration221": 1, "iteration222": 1, "iteration223": 1, "iteration224": 1, "iteration225": 1, "iteration226": 1, "iteration227": 1, "iteration228": 1, "iteration229": 1, "iteration230": 1, "iteration231": 1, "iteration232": 1, "iteration233": 1, "iteration234": 1, "iteration235": 1, "iteration236": 1, "iteration237": 1, "iteration238": 1, "iteration239": 1, "iteration240": 1, "iteration241": 1, "iteration242": 1, "iteration243": 1, "iteration244": 1, "iteration245": 1, "iteration246": 1, "iteration247": 1, "iteration248": 1, "iteration249": 1, "iteration250": 1, "iteration251": 1, "iteration252": 1, "iteration253": 1, "iteration254": 1, "iteration255": 1, "iteration256": 1, "iteration257": 1, "iteration258": 1, "iteration259": 1, "iteration260": 1, "iteration261": 1, "iteration262": 1, "iteration263": 1, "iteration264": 1, "iteration265": 1, "iteration266": 1, "iteration267": 1, "iteration268": 1, "iteration269": 1, "iteration270": 1, "iteration271": 1, "iteration272": 1, "iteration273": 1, "iteration274": 1, "iteration275": 1, "iteration276": 1, "iteration277": 1, "iteration278": 1, "iteration279": 1, "iteration280": 1, "iteration281": 1, "iteration282": 1, "iteration283": 1, "iteration284": 1, "iteration285": 1, "iteration286": 1, "iteration287": 1, "iteration288": 1, "iteration289": 1, "iteration290": 1, "iteration291": 1, "iteration292": 1, "iteration293": 1, "iteration294": 1, "iteration295": 1, "iteration296": 1, "iteration297": 1, "iteration298": 1, "iteration299": 1, "iteration300": 1, "user": 1, "system": 1, "elaps": 1, "604": 1, "020": 1, "120": 1, "926": 1, "30334101718324": 1, "00260494982243339": 1, "00204674530626065": 1, "30096404267195": 1, "00907076937971787": 1, "29863432626151": 1, "0177615086442737": 1, "29307714070162": 1, "0130000283642092": 1, "30428848256149": 1, "00702147549591139": 1, "first": 2, "consid": 2, "simpl": 2, "comput": 2, "extens": 2, "load": [2, 4], "To": 2, "reload": 2, "us": [2, 3, 5], "reload_ext": 2, "warn": [2, 3], "pyplot": [2, 3, 4], "plt": [2, 3, 4], "bilbi": [2, 3, 4], "out_lin": 2, "filterwarn": 2, "ignor": 2, "getlogg": [2, 4], "setlevel": [2, 4], "critic": [2, 4], "injection_paramet": [2, 4], "sampling_frequ": [2, 4], "time_dur": 2, "arang": 2, "ax": [2, 4], "subplot": [2, 4], "o": [2, 5], "label": [2, 3, 4], "signal": [2, 4], "set_xlabel": [2, 4], "y": 2, "class": [2, 3], "gaussianlikelihood": 2, "analytical1dlikelihood": 2, "__init__": [2, 3], "self": [2, 3], "func": 2, "kwarg": [2, 3], "super": [2, 3], "residu": 2, "core": [2, 3], "priordict": [2, 3], "linear_regress": 2, "res_fn": 2, "_result": [2, 4], "json": [2, 4], "read_in_result": 2, "run_sampl": [2, 3, 4], "dynesti": [2, 3, 4], "nlive": [2, 3], "1500": 2, "log_evid": [2, 3, 4], "log_evidence_err": [2, 3, 4], "log_bayes_factor": [2, 4], "08": [2, 3], "nan": [2, 3], "371619": 2, "960707": 2, "144": 2, "432015": 2, "995732": 2, "600449": 2, "023405": 2, "965399": 2, "630136": 2, "349366": 2, "885341": 2, "410890": 2, "399176": 2, "800271": 2, "537964": 2, "386132": 2, "491811": 2, "plot_corn": [2, 4], "suptitl": 2, "lnz_file": [2, 4], "point": [2, 3], "4050": 2, "best_lnz": 2, "best_sampl": 2, "corner": 2, "color": [2, 4], "tab": [2, 4], "grai": 2, "truth": 2, "label_kwarg": 2, "fontsiz": 2, "quantil": 2, "overplot": 2, "enumer": 2, "overplot_lin": 2, "overplot_point": 2, "marker": 2, "markers": 2, "plot_ci": [2, 4], "ravel": [2, 4], "c1": 2, "axhlin": 2, "median": [2, 3], "set_xlim": [2, 4], "set_ylim": [2, 4], "147": 2, "tight_layout": [2, 4], "tmp": 2, "ipykernel_592353": 2, "845767860": 2, "userwarn": 2, "figur": 2, "layout": 2, "ha": 2, "chang": 2, "tight": 2, "145": 2, "hi": 2, "here": 2, "two": 2, "vs": 2, "curv": 2, "sampling_lnz_lnz_err": 2, "figsiz": [2, 4], "xscale": 2, "xlim": 2, "ylim": 2, "150": 2, "xlabel": 2, "ylabel": 2, "ln": [2, 4], "z": 2, "fill_between": 2, "y1": 2, "min": 2, "y2": 2, "interpol": 2, "zorder": 2, "frameon": 2, "3923": 2, "lartillotlikelihood": 3, "dim": 3, "rang": 3, "get_prior": 3, "pri": 3, "truncatednorm": 3, "true_lnz": 3, "nested_sampling_lnz": 3, "lartillot_dynesty_d": 3, "_v": 3, "rwalk": 3, "parallel_tempering_lnz": 3, "bilby_mcmc": 3, "ntemp": 3, "lartillot_ptmc_d": 3, "proposal_cycl": 3, "default": 3, "printdt": 3, "nsampl": 3, "__simulate_posterior": 3, "nsamp": 3, "updat": 3, "ln_prob": 3, "fi_lnz": 3, "linspac": 3, "90": 3, "std": 3, "print_result_dict": 3, "lartillot": 3, "isinst": 3, "tupl": 3, "results_1d": 3, "nest": 3, "parallel": 3, "temper": 3, "06": 3, "results_20d": 3, "run": [3, 4], "lartillot_dynesty_d20_v0": 3, "output": 3, "save": 3, "cach": 3, "pypoetri": 3, "virtualenv": 3, "zrw66etn": 3, "py3": 3, "python3": 3, "site": 3, "util": 3, "73": [3, 4], "deprecationwarn": 3, "access": 3, "attr": 3, "__version__": 3, "deprec": 3, "futur": 3, "releas": 3, "importlib": 3, "metadata": 3, "directli": 3, "queri": 3, "vdict": 3, "getattr": 3, "sy": 3, "modul": 3, "analysi": 3, "x0": 3, "x1": 3, "x2": 3, "x3": 3, "x4": 3, "x5": 3, "x6": 3, "x7": 3, "x8": 3, "x9": 3, "x10": 3, "x11": 3, "x12": 3, "x13": 3, "x14": 3, "x15": 3, "x16": 3, "x17": 3, "x18": 3, "x19": 3, "__main__": 3, "singl": 3, "took": 3, "886e": 3, "bound": 3, "live": 3, "update_interv": 3, "600": 3, "first_upd": 3, "npdim": 3, "rstate": 3, "queue_siz": 3, "pool": 3, "use_pool": 3, "live_point": 3, "logl_arg": 3, "logl_kwarg": 3, "ptform_arg": 3, "ptform_kwarg": 3, "gradient": 3, "grad_arg": 3, "grad_kwarg": 3, "compute_jac": 3, "enlarg": 3, "bootstrap": 3, "walk": 3, "facc": 3, "slice": 3, "fmove": 3, "max_mov": 3, "update_func": 3, "ncdim": 3, "blob": 3, "save_histori": 3, "history_filenam": 3, "maxit": 3, "maxcal": 3, "dlogz": [3, 4], "logl_max": 3, "inf": 3, "n_effect": 3, "add_liv": 3, "print_progress": 3, "print_func": 3, "_print_func": 3, "object": 3, "0x7f9d97721f10": 3, "save_bound": 3, "checkpoint_fil": 3, "checkpoint_everi": 3, "resum": 3, "checkpoint": 3, "everi": 3, "check_point_delta_t": 3, "version": 3, "implement": [3, 5], "act": 3, "averag": 3, "step": [3, 4], "accept": 3, "up": 3, "chain": 3, "length": 3, "5000": 3, "gener": 3, "initi": 3, "written": 3, "01_resum": 3, "pickl": 3, "reject": 3, "obtain": 3, "7823": 3, "59": [3, 4], "191276": 3, "summari": 3, "ln_noise_evid": 3, "256": 3, "223": 3, "lartillot_ptmc_d20_v0": 3, "657e": 3, "nensembl": 3, "pt_ensembl": 3, "tmax": 3, "tmax_from_snr": 3, "initial_beta": 3, "adapt": 3, "adapt_t0": 3, "adapt_nu": 3, "pt_rejection_sampl": 3, "burn_in_nact": 3, "thin_by_nact": 3, "fixed_discard": 3, "autocorr_c": 3, "l1step": 3, "l2step": 3, "1800": 3, "min_tau": 3, "stop_after_converg": 3, "fixed_tau": 3, "tau_window": 3, "evidence_method": 3, "stepping_ston": 3, "initial_sample_method": 3, "initial_sample_dict": 3, "bilbyptmcmcsampl": 3, "converg": 3, "convergenceinput": 3, "target_nsampl": 3, "paralleltemperinginput": 3, "input": 3, "71687116": 3, "51390427": 3, "36840315": 3, "2640976": 3, "18932395": 3, "13572088": 3, "09729439": 3, "06974754": 3, "2795374292990265": 3, "1717320863039971": 3, "05968361770550555": 3, "5230902286641548": 3, "6633760314925232": 3, "21294322516342162": 3, "585494068660412": 3, "8095976482225977": 3, "9182833389630245": 3, "47189605378706473": 3, "4904898766971261": 3, "4137114591843285": 3, "3125255400114014": 3, "3159452862583118": 3, "3666388720970858": 3, "407628354348433": 3, "8250533346649914": 3, "79153659896676": 3, "693913217753955": 3, "0375215193741993": 3, "proposalcycl": 3, "adaptivegaussianpropos": 3, "acceptance_ratio": 3, "differentialevolutionpropos": 3, "uniformpropos": 3, "kdepropos": 3, "train": 3, "fishermatrixpropos": 3, "gmmpropos": 3, "convergence_input": 3, "draw": 3, "80e": 3, "60e": 3, "ad": 3, "e": [3, 5], "0e": 3, "10m": 3, "ev": 3, "maxl": 3, "662": 3, "etf": 3, "50e": 3, "40e": 3, "15m": 3, "39e": 3, "28m": 3, "85e": 3, "84e": 3, "21m": 3, "31e": 3, "29e": 3, "18m": 3, "77e": 3, "19m": 3, "24e": 3, "22e": 3, "70e": 3, "20m": 3, "16e": 3, "15e": 3, "64e": 3, "63e": 3, "09m": 3, "11e": 3, "17m": 3, "57e": 3, "56e": 3, "05e": 3, "04e": 3, "13m": 3, "53e": 3, "52e": 3, "11m": 3, "98e": 3, "97e": 3, "24m": 3, "45e": 3, "16m": 3, "92e": 3, "90e": 3, "14m": 3, "38e": 3, "12m": 3, "88e": 3, "86e": 3, "34e": 3, "81e": 3, "79e": 3, "temperatur": 3, "finish": 3, "03e": 3, "08e": 3, "13e": 3, "08m": 3, "18e": 3, "27e": 3, "32e": 3, "36e": 3, "41e": 3, "46e": 3, "51e": 3, "55e": 3, "65e": 3, "75e": 3, "07m": 3, "89e": 3, "94e": 3, "99e": 3, "07": 3, "03m": 3, "02m": 3, "23e": 3, "06m": 3, "28e": 3, "33e": 3, "05m": 3, "43e": 3, "48e": 3, "62e": 3, "67e": 3, "72e": 3, "76e": 3, "91e": 3, "96e": 3, "01e": 3, "06e": 3, "04m": 3, "10e": 3, "20e": 3, "25e": 3, "30e": 3, "44e": 3, "49e": 3, "54e": 3, "59e": 3, "68e": 3, "82e": 3, "87e": 3, "35e": 3, "sklearn": 3, "1151": [3, 4], "convergencewarn": 3, "number": 3, "distinct": 3, "cluster": 3, "found": [3, 4], "smaller": 3, "n_cluster": 3, "possibl": 3, "due": 3, "duplic": 3, "fit_method": 3, "arg": 3, "29m": 3, "30m": 3, "fail": 3, "refit": 3, "error": 3, "th": 3, "lead": 3, "minor": 3, "posit": 3, "definit": 3, "singular": 3, "01m": 3, "71m": 3, "58e": 3, "66m": 3, "73m": 3, "61e": 3, "results_100d": 3, "usageerror": 4, "magic": 4, "out_gwinj": 4, "bilby_logg": 4, "res_fnam": 4, "durat": 4, "sampling_freq": 4, "min_freq": 4, "1024": 4, "mass": 4, "spin": 4, "ra": 4, "375": 4, "dec": 4, "2108": 4, "extrins": 4, "659": 4, "1126259642": 4, "413": 4, "inj_m1": 4, "inj_m2": 4, "inj_chirp_mass": 4, "convers": 4, "component_masses_to_chirp_mass": 4, "inj_q": 4, "component_masses_to_mass_ratio": 4, "waveform_gener": 4, "waveformgener": 4, "frequency_domain_source_model": 4, "sourc": 4, "lal_binary_black_hol": 4, "parameter_convers": 4, "convert_to_lal_binary_black_hole_paramet": 4, "waveform_argu": 4, "waveform_approxim": 4, "imrphenomd": 4, "reference_frequ": 4, "minimum_frequ": 4, "detector": 4, "ligo": 4, "hanford": 4, "h1": 4, "design": 4, "sensit": 4, "ifo": 4, "interferometerlist": 4, "set_strain_data_from_power_spectral_dens": 4, "start_tim": 4, "inject_sign": 4, "chirp": 4, "howev": 4, "ar": 4, "quit": 4, "un": 4, "astrophys": 4, "process": 4, "convert": 4, "compon": 4, "bbhpriordict": 4, "perform": 4, "check": 4, "doe": 4, "extend": 4, "space": 4, "longer": 4, "validate_prior": 4, "gravitationalwavetransi": 4, "interferomet": 4, "npoint": 4, "conversion_funct": 4, "generate_all_bbh_paramet": 4, "result_class": 4, "cbcresult": 4, "from_json": 4, "filenam": 4, "1984": 4, "70": 4, "store": 4, "so": 4, "instead": 4, "1e1": 4, "2e5": 4, "2653": 4, "bf": 4, "69": 4, "pre": 4, "link": 4, "http": 4, "org": 4, "record": 4, "7884973": 4, "128_cpus_v1_0_result": 4, "hdf5": 4, "from_hdf5": 4, "search_parameter_kei": 4, "15196": 4, "16332": 4, "68": 4, "96249503300896": 4, "303518": 4, "137896": 4, "406246": 4, "221114": 4, "224999": 4, "640732": 4, "363822": 4, "475442": 4, "716670": 4, "535031": 4, "122786": 4, "126647": 4, "208230": 4, "124358": 4, "064088": 4, "294212": 4, "671431": 4, "858286": 4, "436716": 4, "239654": 4, "671935": 4, "089454": 4, "899129": 4, "487101": 4, "332991": 4, "617993": 4, "095211": 4, "673787": 4, "440678": 4, "530625": 4, "277061": 4, "083531": 4, "924704": 4, "632946": 4, "622373": 4, "421349": 4, "331844": 4, "584094": 4, "388201": 4, "657011": 4, "cos_theta_jn": 4, "649265": 4, "750138": 4, "906039": 4, "887552": 4, "956641": 4, "033949": 4, "530611": 4, "617934": 4, "037721": 4, "061500": 4, "838546": 4, "858204": 4, "359604": 4, "671526": 4, "573420": 4, "105585": 4, "906103": 4, "371659": 4, "250562": 4, "480921": 4, "552605": 4, "251147": 4, "141375": 4, "281922": 4, "242988": 4, "000104": 4, "000073": 4, "000216": 4, "000176": 4, "000130": 4, "964813": 4, "408287": 4, "456029": 4, "735976": 4, "985682": 4, "725553": 4, "777387": 4, "630868": 4, "647782": 4, "011267": 4, "835": 4, "104547": 4, "644": 4, "254930": 4, "1117": 4, "045547": 4, "937210": 4, "1433": 4, "822042": 4, "1e9": 4, "histogram_fi_lnz": 4, "rval": 4, "rval_threshold": 4, "idx": 4, "argmin": 4, "abov": 4, "hist": 4, "bin": 4, "histtyp": 4, "green": 4, "axvlin": 4, "ls": 4, "pdf": 4, "1e7": 4, "10000000": 4, "1e6": 4, "out_inj2": 4, "u": 5, "rier": 5, "tegratio": 5, "videnc": 5, "calcu": 5, "l": 5, "ation": 5, "python": 5, "marginalis": 5, "rotiroti": 5, "et": 5, "al": 5, "2018": 5, "follow": 5, "demonstr": 5, "how": 5, "calcul": 5, "three": 5, "scenario": 5, "linear": 5, "regress": 5, "inject": 5, "gw150914": 5}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"ex3": 0, "gw150914": 0, "load": 0, "posterior": [0, 4], "nest": [0, 2, 4], "sampl": [0, 2, 4], "lnz": [0, 2, 3, 4], "comput": [0, 4], "fi": [0, 2, 4], "us": [0, 4], "downsampl": [0, 4], "speed": 0, "full": [0, 4], "few": 0, "refernc": 0, "point": [0, 4], "ex2": [1, 4], "lartillot": 1, "s": 1, "model": 1, "ex1": 2, "linear": 2, "regress": 2, "gener": 2, "some": 2, "data": 2, "calcul": [2, 4], "compar": 2, "valu": 2, "comparison": 3, "1d": 3, "case": 3, "20d": 3, "100d": 3, "plot": 3, "cbc": 4, "gw": 4, "inject": 4, "simpl": 4, "2": 4, "paramet": 4, "inferec": 4, "evid": 4, "15": 4, "d": 4, "infer": 4, "fewer": 4, "refer": 4, "funnel": 5}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9, "sphinx": 56}})