#!/usr/bin/env python3

#  fluxes.py
#  
#  Copyright 2016 Christian Diener <mail[at]cdiener.com>
#  
#  MIT license. See LICENSE for more information.

from cobra.solvers import solver_dict, get_solver_name
from cobra.manipulation.modify import convert_to_irreversible
from cobra.io import read_sbml_model, load_matlab_model
from cobra import Reaction
from urllib.request import urlopen
import zipfile
import gzip
import os
import sys
import shutil

HMRURL = "http://www.metabolicatlas.org/assets/hmr/HMRcollection.xml-8849b3803fdcc7dbd26fea46ef0dcc50.zip"

def download_model(url, basename, dir="models"):
    """Downloads a model and changes its compression to gzip"""
    os.makedirs(dir, exist_ok=True)
    basename = os.path.join(dir, basename.replace(" ", "_"))
    
    if os.path.exists(basename + ".xml.gz"): return
    
    with urlopen(url) as response, open(basename + ".zip", "wb") as f_out:
            shutil.copyfileobj(response, f_out)
    
    with zipfile.ZipFile(basename + ".zip") as zipf:
        xml_files = [z for z in zipf.namelist() if z.endswith(".xml")]
        if len(xml_files) != 1: 
            raise ValueError("Zip file must contain exactly one XML file!")
        with zipf.open(xml_files[0]) as f_in, \
            gzip.open(basename + ".xml.gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    os.remove(basename + ".zip")

def setup_solver(cobra_model, vmax=1.0e16, solver="cglpk", **optimize_kwargs):
    """setup an LP solver instance for pFBA that can be recycled"""
    obj_ids = []
    
    solver = solver_dict[get_solver_name() if solver is None else solver]
    lp = solver.create_problem(cobra_model, **optimize_kwargs)
    
    for i, r in enumerate(cobra_model.reactions):
        solver.change_variable_objective(lp, i, 1.0)
        solver.change_variable_bounds(lp, i, 0.0, vmax)
        if r.objective_coefficient != 0: obj_ids.append(i)
   
    if len(obj_ids) != 1:
       raise ValueError("Need exactly one objective reaction!")
       
    return solver, lp, obj_ids[0]

def to_reversible(irrev_ids, f_dict):
    for rid in irrev_ids:
        f_dict[rid] -= f_dict.pop(rid + "_reverse", 0)
    
    return f_dict

def minimize_total_flux(lp, solver, obj_id, vbm, model, **solver_args):
    """calculate minimum total flux"""
    old = sys.stdout
    f = open(os.devnull, 'w')
    sys.stdout = f
    solver.change_variable_bounds(lp, obj_id, vbm, vbm)
    sol = solver.solve_problem(lp, objective_sense="minimize", **solver_args)
    sys.stdout = old
    f.close()
    if sol != 'optimal':
        raise ValueError("Solution status is: {}!".format(sol))
    
    sol = solver.format_solution(lp, model)
    return sol.x_dict
    
if __name__ == "__main__":
    import pandas as pd
    from timeit import default_timer as timer
    from collections import defaultdict
    
    start = timer()
    
    tissues = pd.read_csv("tissues.csv").set_index("panel")
    pred = pd.read_csv("pred_rates.csv")
    pred = pred[pred["panel"].isin(tissues.index) & pred["tumor"]]
    
    # Download cancer models from HMA
    tissues[["url","tissue"]].drop_duplicates().apply(
        lambda x: download_model(*x), axis=1)
    
    # Calculate pFBA fluxes
    n = pred.shape[0]
    fluxes = []
    samples = []
    subsystem = defaultdict(lambda: "None")
    co = 0
    
    for p in tissues.index.unique():
        path = os.path.join("models", tissues.loc[p, "tissue"] + ".xml.gz")
        model = read_sbml_model(path)
        for r in model.reactions: 
            if subsystem[r.id] != "None" and subsystem[r.id] != r.subsystem:
                raise ValueError("Conflicting subsystem info!")
            subsystem[r.id] = r.subsystem
        irrev_ids = [r.id for r in model.reactions]
        convert_to_irreversible(model)
        model.reactions.get_by_id("CancerBiomass_OF").objective_coefficient = 1
        solver, lp, oid = setup_solver(model)
        data = pred[pred["panel"] == p]
        for i in range(data.shape[0]):
            if co % (n//100) == 0:
                print("                             \r", end="")
                print("Calculated {}% of fluxes...".format(co//(n//100)), end="")
            samples.append(data["barcode"].iloc[i])
            rate = max(data["rates"].iloc[i], 0.0)
            f = minimize_total_flux(lp, solver, oid, rate, model)
            fluxes.append(to_reversible(irrev_ids, f))
            co += 1
    
    print("\nNeeded {:.2f} s.".format(timer() - start))
    
    # Save output and subsystems
    print("Saving non-zero fluxes...")
    fluxes = pd.DataFrame(fluxes, index=samples)
    fluxes = fluxes.loc[:, fluxes.abs().max() > 1e-6]
    fluxes.to_csv("fluxes.csv")
    
    print("Saving reaction info...") 
    download_model(HMRURL, "hmr")
    model = read_sbml_model(os.path.join("models", "hmr.xml.gz"))
    rids = fluxes.columns.values
    subsystem = {"id": list(subsystem.keys()), "subsystem": list(subsystem.values())}
    info = pd.DataFrame(subsystem)
    
    info.to_csv("flux_info.csv")
        
