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
import shutil

HMRURL = "http://www.metabolicatlas.org/assets/hmr/HMRcollection2_00.xml-f0de1f951d16f78abf131cece19f8af7.zip"

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

#def build_objective(model, ref):
    #rbm = ref.reactions.get_by_id("biomass_components")
    #mets = [[m for m in model.metabolites if m.compartment == k.compartment and \
            #m.name == k.name] for k in rbm.metabolites.keys()]
    #v = list(rbm.metabolites.values())
    #nm = [len(m) for m in mets]
    #if max(nm) == 0: raise ValueError("Could not map reactions!")
    #if max(nm) > 1: raise ValueError("Duplicate metabolite names in mapping!")
    #stoich = { mets[i][0]: v[i] for i in range(len(v)) if nm[i] > 0 }
    #rnew = Reaction("biomass_components")
    #rnew.add_metabolites(stoich)
    #rnew.lower_bound = 0
    #rnew.upper_bound = 1.0e16
    #rnew.objective_coefficient = 1.0
    #return rnew

def to_reversible(irrev_ids, f_dict):
    for rid in irrev_ids:
        f_dict[rid] -= f_dict.pop(rid + "_reverse", 0)
    
    return [f_dict[k] for k in irrev_ids]

def minimize_total_flux(lp, solver, obj_id, vbm, model, **solver_args):
    """calculate minimum total flux"""
    solver.change_variable_bounds(lp, obj_id, vbm, vbm)
    sol = solver.solve_problem(lp, objective_sense="minimize", **solver_args)
    if sol != 'optimal':
        raise ValueError("Solution status is: {}!".format(sol))
    
    sol = solver.format_solution(lp, model)
    return sol.x_dict
    
if __name__ == "__main__":
    import pandas as pd
    from timeit import default_timer as timer
    import sys
    
    start = timer()
    
    model_file = sys.argv[1]
    _, ext = os.path.splitext(model_file)
    
    pred = pd.read_csv("pred_rates.csv")
    
    if ext.startswith(".xml"):
        model = read_sbml_model(model_file)
    elif ext == ".mat":
        model = load_matlab_model(model_file)
    irrev_ids = [r.id for r in model.reactions]
    convert_to_irreversible(model)
    
    # Download tissue models
    n = pred.shape[0]
    fluxes = []
    
    solver, lp, oid = setup_solver(model)
    for i in range(n):
        if i % (n//100) == 0:
            print("                             \r", end="")
            print("Calculated {}% of fluxes...".format(i//(n//100)), end="")
        rate = max(pred["rates"].iloc[i], 0.0)
        f = minimize_total_flux(lp, solver, oid, rate, model)
        fluxes.append(to_reversible(irrev_ids, f))
    print("\nNeeded {:.2f} s.".format(timer() - start))
    print("Saving non-zero fluxes...")
    fluxes = pd.DataFrame(fluxes, columns=irrev_ids)
    fluxes = fluxes.loc[:, fluxes.abs().max() > 1e-6]
    fluxes.to_csv("fluxes.csv") 
    rids = fluxes.columns.values
    info = pd.DataFrame({"id": rids, 
        "name": [model.reactions.get_by_id(i).name for i in rids],
        "reaction": [model.reactions.get_by_id(i).reaction for i in rids],
        "subsystem": [model.reactions.get_by_id(i).subsystem for i in rids]})
    info.to_csv("flux_info.csv")
        
