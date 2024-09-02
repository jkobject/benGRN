import os
from pathlib import Path
from typing import Dict

# Local imports
from . import run_experiment as runexp
from . import post_processing as po


class BoolODE(object):
    """
    The BoolODE object is created by parsing a user-provided configuration
    file. It then calls the correct functions for executing the user defined
    experiment.
    """

    def __init__(
        self,
        jobs,  # joblist {"name": "dyn-LL-1", "simulation_time": 9, "num_cells": 300, "do_parallel": True, "sample_cells": False, 'nClusters': 1}
        dropout_jobs: Dict[
            str, str
        ] = None,  # Dropouts {"drop_cutoff": 0.5, "drop_prob": 0.4, "sample_size": 100}
        perplexities=None,  #
        slingshot_jobs=None,  # Slingshot {"perplexity": 200}
        gensample_jobs=None,  # GenSamples {"sample_size":100, "nDatasets":1}
        model_dir="data",
        output_dir="output",
        do_simulations=True,
        do_post_processing=True,
        modeltype="hill",
    ) -> None:
        """
            __init__ _summary_

            Args:
                jobs (_type_): _description_
                output_dir (_type_): _description_
                do_simulations (_type_): _description_
                do_post_processing (_type_): _description_
                modeltype (_type_): _description_
                nCells
        drop_cutoff
        drop_prob
        """
        self.dropout_jobs = dropout_jobs
        self.perplexities = perplexities
        self.slingshot_jobs = slingshot_jobs
        self.gensample_jobs = gensample_jobs
        self.geneexpression_jobs = None
        self.model_dir = model_dir
        self.output_dir = output_dir
        self.do_simulations = do_simulations
        self.do_post_processing = do_post_processing
        self.modeltype = modeltype
        self.jobs = self._process_jobs([jobs])

    def _process_jobs(self, jobs) -> Dict[int, Dict]:
        """
        Creates a list of jobs, where each job is specified by
        the values in a dictionary called data.
        Default parameter values are specified here.
        """
        n_jobs = {}
        for jobid, job in enumerate(jobs):
            data = {}
            # Create output folder if it doesnt exist
            data["name"] = job.get("name")
            data["outprefix"] = Path(self.output_dir, job.get("name"))
            data["modelpath"] = Path(self.model_dir, job.get("model_definition"))
            data["simulation_time"] = job.get("simulation_time", 20)
            data["icsPath"] = Path(
                self.model_dir, job.get("model_initial_conditions", "")
            )
            data["num_cells"] = job.get("num_cells", 100)
            data["sample_cells"] = job.get("sample_cells", False)
            data["nClusters"] = job.get("nClusters", 1)
            data["doParallel"] = job.get("do_parallel", False)
            data["identical_pars"] = job.get("identical_pars", False)
            data["sample_pars"] = job.get("sample_pars", False)
            data["sample_std"] = job.get("sample_std", 0.1)
            data["integration_step_size"] = job.get("integration_step_size", 0.01)
            # Optional Settings
            data["parameter_inputs_path"] = Path(
                self.model_dir, job.get("parameter_inputs_path", "")
            )
            data["parameter_set"] = Path(self.model_dir, job.get("parameter_set", ""))
            data["interaction_strengths"] = Path(
                self.model_dir, job.get("interaction_strengths", "")
            )
            data["species_type"] = Path(self.model_dir, job.get("species_type", ""))
            # Simulator settings
            data["burnin"] = job.get("burn_in", False)
            data["writeProtein"] = job.get("write_protein", False)
            data["normalizeTrajectory"] = job.get("normalize_trajectory", False)
            data["add_dummy"] = job.get("add_dummy", False)
            data["max_parents"] = job.get("max_parents", 1)
            data["modeltype"] = self.modeltype

            n_jobs[jobid] = data
        return n_jobs

    def execute_jobs(self, parallel=False, num_threads=1):
        """
        Run each user specified job.
        BoolODE runs two types of functions
        1. If `do_simulation == TRUE`, perform SDE simulations of model specified as Boolean rules.
        2. If `do_post_processing == TRUE` perform the list of post processing operations specified.

        .. warning::
            This function automatically creates folders for each job name
            as specified in the config file, if the folder doesn't already exist.
            Contents of existing folders will be rewritten!
        """
        _ = self.output_dir

        alljobs = self.jobs.keys()
        print("Creating output folders")
        for jobid in alljobs:
            outdir = self.jobs[jobid]["outprefix"]
            if not os.path.exists(outdir):
                print(outdir, "does not exist, creating it...")
                os.makedirs(outdir)
        if self.do_simulations:
            print("Starting simulations")
            for jobid in alljobs:
                runexp.startRun(self.jobs[jobid])
        if self.do_post_processing:
            print("Starting post processing")
            self.do_post_processing()

    def do_post_processing(self):
        """
        Call genSamples() first. Then run DimRed runSlingShot,
        generateDropouts, if specified by the user.
        """
        alljobs = list(self.jobs.keys())
        filetypedict = {
            "expr": "ExpressionData.csv",
            "pseudo": "PseudoTime.csv",
            "refNet": "refNetwork.csv",
        }

        ## Always do genSamples() once if even a single other analysis is requested!
        doOtherAnalysis = False
        if (
            self.dropout_jobs
            or self.perplexities
            or self.geneexpression_jobs
            or self.slingshot_jobs
        ):
            doOtherAnalysis = True
        generatedPaths = {}

        if self.gensample_jobs is not None or doOtherAnalysis:
            print("Generating Samples...")
            if self.gensample_jobs is None:
                gsamp = {}
                gsamp["sample_size"] = self.jobs[alljobs[0]]["num_cells"]
                gsamp["nDatasets"] = 1
                gensample_jobs = [gsamp]
            else:
                gensample_jobs = self.gensample_jobs
            for gsamp in gensample_jobs:
                for jobid in alljobs:
                    settings = {}
                    settings["num_cells"] = self.jobs[jobid]["num_cells"]
                    settings["sample_size"] = gsamp.get("sample_size", 100)
                    settings["outPrefix"] = str(self.jobs[jobid]["outprefix"])
                    settings["nDatasets"] = gsamp.get("nDatasets", 1)
                    settings["name"] = self.jobs[jobid]["name"]
                    settings["nClusters"] = self.jobs[jobid]["nClusters"]
                    generatedPaths[jobid] = po.genSamples(settings)

        if self.dropout_jobs is not None:
            print("Starting genDropouts...")
            for drop in self.dropout_jobs:
                num_invalid = 0
                for jobid in alljobs:
                    for gsampPath in generatedPaths[jobid]:
                        settings = {}
                        invalid = False
                        settings["outPrefix"] = gsampPath
                        settings["expr"] = Path(gsampPath, "ExpressionData.csv")
                        settings["pseudo"] = Path(gsampPath, "PseudoTime.csv")
                        settings["refNet"] = Path(gsampPath, "refNetwork.csv")
                        settings["dropout"] = drop.get("dropout", True)
                        settings["sample_size"] = drop.get("sample_size", 100)
                        settings["num_cells"] = self.jobs[jobid]["num_cells"]
                        settings["drop_cutoff"] = drop.get("drop_cutoff", 0.0)
                        settings["drop_prob"] = drop.get("drop_prob", 0.0)

                        for filetype in ["expr", "pseudo", "refNet"]:
                            if not settings[filetype].is_file():
                                print(
                                    self.jobs[jobid]["name"],
                                    ": ",
                                    filetypedict[filetype],
                                    "not found. Retry with `do_simulations: True` in global_settings.",
                                )
                                invalid = True
                                num_invalid += 1
                                break
                        if not invalid:
                            po.genDropouts(settings)
                    if num_invalid == len(alljobs):
                        break

        if self.perplexities is not None:
            print("Starting dimesionality reduction using tSNE")
            for perplexity in self.perplexities:
                num_invalid = 0
                for jobid in alljobs:
                    for gsampPath in generatedPaths[jobid]:
                        print("perplexity=", perplexity)
                        settings = {}
                        invalid = False
                        settings["expr"] = Path(gsampPath, "ExpressionData.csv")
                        settings["pseudo"] = Path(gsampPath, "PseudoTime.csv")
                        settings["perplexity"] = perplexity
                        settings["default"] = False
                        for filetype in ["expr", "pseudo"]:
                            if not settings[filetype].is_file():
                                print(
                                    self.jobs[jobid]["name"],
                                    ":",
                                    filetypedict[filetype],
                                    "not found. Retry with `do_simulations: True` in global_settings.",
                                )
                                invalid = True
                                num_invalid += 1
                                break
                        if not invalid:
                            po.doDimRed(settings)
                            _ = False
                    if num_invalid == len(alljobs):
                        break

        if self.geneexpression_jobs is not None:
            if self.perplexities is None:
                print("Using default perplexity=50 (Specify `perplexity` under DimRed)")
                perplexity = 50
            else:
                if len(self.perplexities) > 1:
                    perplexity = min([j["perplexity"] for j in self.perplexities])
                else:
                    perplexity = self.perplexities[0]["perplexity"]

            print("Plotting gene expression levels in tSNE projection")
            for jobid in alljobs:
                for gsampPath in generatedPaths[jobid]:
                    print(gsampPath)
                    settings = {}
                    invalid = False
                    settings["expr"] = Path(gsampPath, "ExpressionData.csv")
                    settings["pseudo"] = Path(gsampPath, "PseudoTime.csv")
                    settings["perplexity"] = perplexity
                    settings["default"] = False
                    po.plotGeneExpression(settings)

        if self.slingshot_jobs is not None:
            if self.perplexities is None:
                print(
                    "Using default perplexity=300. (Specify `perplexity` under DimRed.)"
                )
            print("Starting SlingShot...")
            for sshot in self.slingshot_jobs:
                for jobid in alljobs:
                    for gsampPath in generatedPaths[jobid]:
                        settings = {}
                        invalid = False
                        settings["outPrefix"] = (
                            gsampPath + "/" + gsampPath.split("/")[-1] + "-ss"
                        )
                        settings["expr"] = Path(gsampPath, "ExpressionData.csv")
                        settings["pseudo"] = Path(gsampPath, "PseudoTime.csv")
                        settings["refNet"] = Path(gsampPath, "refNetwork.csv")
                        if self.jobs[jobid]["nClusters"] == 1:
                            settings["nClusters"] = 1
                        else:
                            settings["nClusters"] = self.jobs[jobid]["nClusters"] + 1
                        settings["noEnd"] = sshot.get("noEnd", False)
                        settings["perplexity"] = sshot.get("perplexity", 300)

                        for filetype in ["expr", "pseudo", "refNet"]:
                            if not settings[filetype].is_file():
                                print(
                                    self.jobs[jobid]["name"],
                                    ": ",
                                    filetypedict[filetype],
                                    "not found. Retry with `do_simulations: True` in global_settings.",
                                )
                                invalid = True
                                break
                        if not invalid:
                            po.computeSSPT(settings)
