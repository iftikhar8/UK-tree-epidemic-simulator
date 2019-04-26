parameters = {"time": 10, "time_horizon": 5000, "ensemble": range(ensemble_size), "t_init": [5, 6], "L": 500}
settings = {"out_path": output_path, "domain_type": domain_type, "date": date, "job_id": job_id, "param_dim": 10,
            "ensemble": range(ensemble_size), "metrics": ["eff", "chi", "d", "t", "P"], "plt_tseries": False,
            "save_figs": False, "dyn_plts": [False, 1, True], "anim": False, "BCD3": True,
            "felling": [False, "40_50", 0.50]}

