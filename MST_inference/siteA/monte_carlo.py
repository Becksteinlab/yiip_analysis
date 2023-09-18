#!/usr/bin/env python

from multibind.multibind import MultibindScanner
from pathlib import Path
import numpy as np
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from collections import defaultdict
import random

# Script to run the finite temperature monte carlo for the thermodyanamical model.

class Stepper(object):

    def __init__(self, groups, rundir, deventries=20, dev_max=0.1, dev_min=1e-3):

        self.logfilename = rundir / "steps.log"
        self.deventries = deventries

        # is this a continued run?
        if self.logfilename.exists():
            self.read_log()
        else:
            self.deviations = []
            self.success_count = 0
            self.total_count = 0

        self.logger = logging.getLogger('MC')
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler(rundir / "steps.log")
        self.logger.addHandler(fh)
        self.load_MST_data()
        self.groups = groups
        self.imgdir = rundir / "img"
        self.inputdir = rundir / "inputs"
        self.repeated_failures = 0
        self.dev_max = dev_max
        self.dev_min = dev_min

        input_states = self.inputdir / "states.csv"
        self.input_graph_best = self.inputdir / "best_graph.csv"

        if not self.input_graph_best.exists():
            input_graph = self.inputdir / "graph.csv"
        else:
            input_graph = self.input_graph_best

        self.scanner = MultibindScanner(input_states, input_graph)

    def read_log(self):
        data = []
        with open(self.logfilename, 'r') as F:
            for line in F:
                if line.startswith("STEP"):
                    data.append(line.strip())
        last_entry = data[-1].split(" ")
        self.total_count = int(last_entry[0].split("=")[1])
        self.success_count = int(last_entry[1].split("=")[1])
        self.deviations = [float(last_entry[5].split("=")[1])] * self.deventries

    def proton_uptake(self):

        msp = self.scanner.results.microstate_probs

        pH = [5.6, 6.0, 6.5, 7.0, 7.4]

        expected = defaultdict()

        for p in pH:
            prob = np.dot(np.array(self.scanner.c.states.Zn), msp.sel({'pH': p}).values)
            expected[p] = prob

        return self.scanner.results['Zn'].values, expected

    def load_MST_data(self):
        mst_data_dir = Path("target")
        pHs = [5.6, 6.0, 6.5, 7.0, 7.4]
        data = defaultdict()
        for ph in pHs:
            df = pd.read_csv(mst_data_dir / 'pH{}.csv'.format(ph))
            data[ph] = df

        self.target = data

    def create_pH_range(self):

        data = self.target

        pH = [5.6, 6.0, 6.5, 7.0, 7.4]
        self.pH = pH

        concentrations = [np.array([1])]
        for ph in pH:
            df = data[ph]
            concentrations.append(df.Dose)
        
        concentrations = np.concatenate(np.array(concentrations))
        concentrations = np.unique(concentrations)
        concentrations.sort()

        masks = defaultdict(list)

        for concentration in concentrations:
            for ph in pH:
                if concentration in data[ph].Dose.values:
                    masks[ph].append(True)
                else:
                    masks[ph].append(False)


        self.masks = masks
        self.concentrations = {'pH': self.pH, 'Zn': concentrations}

    def results_rmsd(self):

        masks = self.masks
        data = self.target

        scanner_zn, uptake = self.proton_uptake()

        mst_values = []
        filtered_uptake = []
        for p in self.pH:
            mst_values.append(data[p].Response.values)
            filtered_uptake.append(uptake[p][masks[p]])
        mst_values = np.concatenate(mst_values)
        filtered_uptake = np.concatenate(filtered_uptake)

        RMSD = np.sqrt(np.mean((mst_values - filtered_uptake)**2))

        return RMSD

    def step(self):

        scanner = self.scanner

        old_graph = scanner.c.graph.copy()
        old_rmsd = self.results_rmsd()
        old_results = scanner.results.copy()

        standard_deviation = self.dev_max

        deviations = np.random.uniform(-0.2, 0.2, size=len(self.groups))
        for deviation, group in zip(deviations, self.groups):

            scanner.c.graph.loc[group[0][2] - 1, ['value']] += deviation

            for g in group:
                scanner.c.graph.loc[g[2] - 1, ['value']] = scanner.c.graph.loc[group[0][2] - 1, ['value']]

        scanner.run(self.concentrations)

        new_rmsd = self.results_rmsd()
        rmsd = self.results_rmsd()

        success = False
        drmsd = new_rmsd-old_rmsd
        p = random.random()
        p_ref = np.exp(-drmsd/0.0001)

        if p <= p_ref:
            for i in scanner.c.graph.iterrows():
                s1 = i[1].state1
                s2 = i[1].state2
                if ('Zn' in s1) ^ ('Zn' in s2):
                    dG = float(scanner.results.dGs.sel(pH=7.0, Zn=1.0, state_i=s1, state_j=s2).values)

                    scanner.c.graph.loc[i[0], "value"] = dG
                    scanner.c.graph.loc[i[0], "variance"] = 1.0
                else:
                    dG = float(scanner.results.dGs.sel(pH=7.0, Zn=1.0, state_i=s1, state_j=s2).values)
                    pKa = 7 - dG / np.log(10)

                    scanner.c.graph.loc[i[0], "value"] = pKa
                    scanner.c.graph.loc[i[0], "variance"] = 1.0

            success = True
            self.success_count += 1
            self.repeated_failures = 0
            self.deviations.append(abs(deviation))
        else:
            scanner.c.graph = old_graph.copy()
            scanner.results = old_results.copy()
            success = False
            self.repeated_failures += 1
        mean_dev = np.mean(self.deviations[-self.deventries:])
        self.logger.info(f"STEP={self.total_count} ACCEPTED={self.success_count} RTOT={rmsd} DRTOT={drmsd} MDEV={mean_dev} STD={standard_deviation}")

        return success, rmsd

    def run(self, target=None, eps=1e-3, allowed_failures=1000, rmsd_ref=0.0765):

        self.create_pH_range()
        self.scanner.run(self.concentrations)

        plt.ion()
        figure, ax = plt.subplots(figsize=[4, 3])
        sns.despine(ax=ax, offset=5)

        plot_ph = defaultdict()
        colors = ['blue', 'orange', 'red', 'green', 'black']
        phs, uptake = self.proton_uptake()
        for ph, color in zip(self.pH, colors):
            uptake_line, = ax.plot(phs[:-1], uptake[ph][:-1], color=color, label=str(ph))
            plot_ph[ph] = uptake_line
            ax.scatter(self.target[ph].Dose, self.target[ph].Response, color=color)

        ax.set_ylabel(r"$\langle x \rangle$")
        ax.set_xlabel(r"Dose")
        ax.set_xscale("log", nonpositive='clip')
        ax.set_ylim([-0.2, 1.2])
        ax.set_xlim([1e-12, 0.05])

        plt.tight_layout()
        plt.savefig(self.imgdir / f"{self.success_count:04d}.png")

        while True:
            success, rmsd = self.step()
            self.total_count += 1

            if not success:
                if self.repeated_failures >= allowed_failures:
                    self.logger.info(f"Failed to find a valid move after {allowed_failures}")
                    break
                continue

            self.scanner.c.write_graph(self.input_graph_best)

            phs, uptake = self.proton_uptake()
            for ph, color in zip(self.pH, colors):
                plot_ph[ph].set_xdata(phs[:-1])
                plot_ph[ph].set_ydata(uptake[ph][:-1])

            figure.canvas.draw()
            figure.canvas.flush_events()

            plt.savefig(self.imgdir / f"{self.success_count:05d}.png")

            if np.mean(self.deviations[-self.deventries:]) < eps and len(self.deviations) >= self.deventries:
                self.logger.info("Converged.")
                break
                
            if rmsd <= rmsd_ref:
                self.logger.info("Converged.")
                break

def main():

    import uuid
    import shutil
    from sys import argv

    runid = uuid.uuid4()

    if len(argv) == 2:
        site = argv[1]

        if site not in ('side', 'center'):
            rundir = Path(site)
            site = str(rundir).split('-')[-1]
        else:
            rundir = Path("runs") / f"{runid}-{site}"

    else:
        site = None
        rundir = Path("runs") / f"{runid}-{site}"

    runinputs = rundir / "inputs"
    runimg = rundir / "img"

    runimg.mkdir(exist_ok=True, parents=True)
    runinputs.mkdir(exist_ok=True, parents=True)

    statefile = Path() / "inputs" / "states.csv"
    graphfile = Path() / "inputs" / "graph.csv"

    if not (runinputs / "states.csv").exists():
        shutil.copy(statefile, runinputs)
    if not (runinputs / "graph.csv").exists():
        shutil.copy(graphfile, runinputs)


    graph = pd.read_csv('inputs/graph.csv')
    groups = []
    for i in graph.iterrows():
        groups.append([(i[1].state1, i[1].state2, i[0]+1)])

    stepper = Stepper(groups, rundir)
    stepper.run(target=site)


if __name__ == "__main__":
    main()