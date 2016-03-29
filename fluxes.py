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
    """Downloads a model and changes its compression to gzip
    
    Downloads a model in zip format, finds a single XML file 
    and converts it to a gzip version. Basically this is done
    since cobrapy can read directly from gzip but not from zip.
    
    Args:
        url: A valid URL string pointing to the location of the
            model.
        basename: A string specifying the new name that will be
            given to the model. Spaces are replaced by underscores.
        dir: The directory where the gzipped model will be saved.
    
    Returns:
        Nothing.
    
    Raises:
        ValueError: The zip file did not contain a XML file. 
    """
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
    """Setup an LP solver instance for pFBA that can be recycled.
    
    Does the initial setup for the pFBA model. Thus, each reaction 
    is set as part of the objective and the bounds are set to [0, Inf]
    where Inf is not really infinite (cobrapy does not allow that) but
    large.
    
    Args:
        cobra_model: A cobra model (class Model) to be set up. Must be
            irreversible. 
        vmax: A float setting the upper bound for reactions. If you want an
            unconstrained model this should be very large.
        solver: The name of the solver to use. 
        optimize_kwargs: Optional solver key word arguments.
    
    Returns:
        A tuple (solver, lp, oid) where solver is the set up solver
        instance, lp the set up LP problem and oid the index of the
        original objective reaction.
        
    Raises:
        ValueError: There was more than one objective reaction.
    """
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

def to_reversible(rev_ids, f_dict):
    """Combines irreversible fluxes to reversible ones.
    
    Args:
        rev_ids: A list of srings denoting the IDs of the reversible
            reactions.
        f_dict: A dictionary where the keys are the irreversible reaction
            IDs and the values the corresponding fluxes.
    
    Returns:
        A dictionary where keys denote reversible reactions IDs and values
        the reversible fluxes.
    """
    for rid in irrev_ids:
        f_dict[rid] -= f_dict.pop(rid + "_reverse", 0)
    
    return f_dict

def minimize_total_flux(solver, lp, obj_id, vbm, model, **solver_args):
    """Calculate minimum total flux.
    
    Obtains the minimum total flux min(sum_i |v_i|) for a given objective value.
    
    Args:
        solver: The solver instance to recycle.
        lp: The LP problem to modify ine each run.
        obj_id: The index of the objective reaction.
        vbm: The upper flux bound.
        model: The model on which to run.
        solver_args: Additional keyword arguments passsed to the solver.
    
    Returns:
        A dictionary {id: val, ...} where id denotes the IDs of irreversible
        reactions and val the flux for that reaction.
        
    Raises:
        ValueError: There is no feasible solution.
    """
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
    
    tissues = pd.read_csv("tissues.csv")
    pred = pd.read_csv("pred_rates.csv")
    pred = pred[pred["panel"].isin(tissues["panel"]) & pred["tumor"]]
    
    # Download cancer models from HMA
    tissues[["url","tissue"]].drop_duplicates().apply(
        lambda x: download_model(*x), axis=1)
    
    # Calculate pFBA fluxes
    n = pred.shape[0]
    fluxes = []
    samples = []
    subsystem = defaultdict(lambda: "None")
    co = 0
    tissues = tissues.set_index("tissue")
    
    for ti in tissues.index.unique():
        path = os.path.join("models", ti + ".xml.gz")
        model = read_sbml_model(path)
        convert_to_irreversible(model)
        for r in model.reactions: 
            if subsystem[r.id] != "None" and subsystem[r.id] != r.subsystem:
                raise ValueError("Conflicting subsystem info!")
            subsystem[r.id] = r.subsystem
        model.reactions.get_by_id("CancerBiomass_OF").objective_coefficient = 1
        solver, lp, oid = setup_solver(model)
        panels = tissues.loc[ti, "panel"]
        if type(panels) == str: panels = [panels] 
        else: panels = panels.tolist() 
        data = pred[pred["panel"].isin(panels)]
        
        for i in range(data.shape[0]):
            if co % (n//100) == 0:
                print("                             \r", end="")
                print("Calculated {}% of fluxes...".format(co//(n//100)), end="")
            co += 1
            rate = data["rates"].iloc[i]
            if rate <= 0: continue
            f = minimize_total_flux(solver, lp, oid, rate, model)
            samples.append(data["barcode"].iloc[i])
            fluxes.append(f)
    
    print("\nNeeded {:.2f} s.".format(timer() - start))
    
    # Save output and subsystems
    print("Saving non-zero fluxes...")
    fluxes = pd.DataFrame(fluxes, index=samples).fillna(0)
    fluxes = fluxes.loc[:, fluxes.abs().max() > 1e-6]
    fluxes.to_csv("fluxes.csv")
    
    print("Saving reaction info...") 
    download_model(HMRURL, "hmr")
    model = read_sbml_model(os.path.join("models", "hmr.xml.gz"))
    rids = fluxes.columns.values
    subsystem = {"id": list(subsystem.keys()), "subsystem": list(subsystem.values())}
    info = pd.DataFrame(subsystem).set_index("id")
    info = info.loc[fluxes.columns]
    
    info.to_csv("flux_info.csv")
        
