# -*- coding: utf-8 -*-
"""
Created on Tue May 30 14:48:03 2023

@author: jonat
"""

# Standard imports.
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import re
import scipy.io
import sys

# Local/third party imports.
from pybasicbayes.util.text import progprint_xrange
from pyhsmm.basic.distributions import PoissonDuration

if not [x for x in sys.path if "pyhsmm_spiketrains" in x]:
    # sys.path.append("/home/jmc/projects/py-hmm/pyhsmm_spiketrains")
    sys.path.append("/home/jmc/Dropbox (NYU Langone Health)/py-hmm-2/pyhsmm_spiketrains")
import pyhsmm_spiketrains.models

reload(pyhsmm_spiketrains.models)
from pyhsmm_spiketrains.internals.match_labels import get_optimal_assignment, get_overlap_mat

# Set whether or not results of run will be saved.
SAVE = True

# Set run date and run ID. Run ID will be used in filename, whereas run date (YYYY/MM/DD) will be used to change to
# save directory (this is mainly to help guard against accidental overwriting). Data directory contains the preprocessed
# data.
run_date = "2023-06-14-A"
run_id = "U"
data_type = "emp/resampled" # "sim" or "emp/resampled"
save_dir = "/home/jmc/Dropbox (NYU Langone Health)/py-hmm-2/results/{}".format(run_date)
data_dir = "/home/jmc/Dropbox (NYU Langone Health)/bhv-mfp/preprocessing/results/preprocessed-data/poisson-hmm/{}".format(data_type)
assert os.path.isdir(save_dir), "Save directory does not exist. You may need to create one with the run_date: {}".format(run_date)


# Set list of subject/sessions as well as training and testing sets for each.
# GC-date-binsize-session part-subsection (subsection should be 0)
# interaction
subj_sess = [[1, '190912', 8, "2ndF1", 0],
             [1, '190912', 8, "2ndM2", 0],
             [6, '200311', 8, "2ndF1", 0], # No usable/worthwhile subsections
             [6, '200311', 8, "2ndM2", 0],
             [6, '200312', 8, "2ndM1", 0],
             [6, '200312', 8, "2ndF2", 0],
             [6, '200312', 8, "2ndT3", 0],
             [8, '200527', 8, "2ndM1", 0],
             [8, '200527', 8, "2ndF2", 0],
             [8, '200611', 8, "2ndF1", 0],
             [8, '200611', 8, "2ndM2", 0],
             [8, '200615', 8, "2ndM1", 0],
             [8, '200615', 8, "2ndF2", 0],
             [8, '200615', 8, "2ndM3", 0],
             [8, '200617', 8, "2ndM1", 0],
             [8, '200617', 8, "2ndF2", 0],
             [12, '200801', 8, "2ndF1", 0],
             [12, '200801', 8, "2ndM2", 0],
             [12, '200801', 8, "2ndT3", 0],
             [12, '200803', 8, "2ndM1", 0],
             [12, '200803', 8, "2ndF2", 0],
             [12, '200803', 8, "2ndT3", 0],
             [17, '200708', 8, "2ndF1", 0],
             [17, '200708', 8, "2ndM2", 0],
             [17, '200709', 7, "2ndF1", 0],
             [17, '200709', 7, "2ndM2", 0],
             [17, '200719', 8, "2ndM1", 0],
             [17, '200719', 8, "2ndF2", 0],
             [17, '200719', 8, "2ndT3", 0],
             [7,  '200527', 8, "2ndM1", 0],
             [7,  '200527', 6, "2ndF2", 0],
             [7,  '200530', 7, "2ndM1", 0],
             [7,  '200530', 7, "2ndF2", 0],
             [7,  '200604', 8, "2ndM1", 0],
             [7,  '200617', 8, "2ndM1", 0],
             [7,  '200617', 8, "2ndF2", 0]]
# wholesession
# subj_sess = [[1, '190912', 8, "wholesession", 0],
#              [6, '200311', 8, "wholesession", 0],
#              [6, '200312', 8, "wholesession", 0],
#              [8, '200527', 8, "wholesession", 0],
#              [8, '200611', 8, "wholesession", 0],
#              [8, '200615', 8, "wholesession", 0],
#              [8, '200617', 8, "wholesession", 0],
#              [12, '200801', 8, "wholesession", 0],
#              [12, '200803', 8, "wholesession", 0],
#              [17, '200708', 8, "wholesession", 0],
#              [17, '200709', 7, "wholesession", 0],
#              [17, '200719', 8, "wholesession", 0],
#              [7,  '200527', 8, "wholesession", 0],
#              [7,  '200530', 7, "wholesession", 0],
#              [7,  '200604', 8, "wholesession", 0],
#              [7,  '200617', 8, "wholesession", 0]]

for zscore_flag in [0]:
    # ------------------ Set loading and model parameters ----------------------- #
    # First, set preprocessing parameters to retrieve desired file.
    pparams = dict(
        data_type="fmpp",
        zscore=zscore_flag,
        relabeled=1,
        resampled=5,
        preprocess_version="B")

    # ----------------- Load and train on training datasets --------------------- #
    # Switch to data directory.
    curr_dir = os.getcwd()


    for i_dataset in range(len(subj_sess)):
        if subj_sess[i_dataset][4] == -1:
            continue

        os.chdir(data_dir)
        load_filename \
            = "{}-{}-{data_type}-zscore{zscore}-bin{}-relabeled{relabeled}-resampled{resampled}-{}-subsection{}-{preprocess_version}".format(
                subj_sess[i_dataset][0],
                subj_sess[i_dataset][1],
                subj_sess[i_dataset][2],
                subj_sess[i_dataset][3],
                subj_sess[i_dataset][4],
                **pparams)
        data = scipy.io.loadmat(load_filename + ".mat")
        neural_data = data["neural_data_binned"]
        neural_data = np.int64(neural_data).transpose()
        labels = data["bin_labels"]

        # Set number of states. This is the Number of states for Dirichlet prior
        # models, or max states for DP prior models.
        n_states = 20

        # Set up model parameters (duration parameters for HSMM models set separately below).
        mparams = dict(
            seed=0,
            model_type="PoissonHDPHSMM",
            n_iter=1000,
            model_args=dict(
                N=np.int(neural_data.shape[1]),  # Number of regions
                alpha_obs=1.0,
                beta_obs=1.0,
                alpha_a_0=10.0,
                alpha_b_0=1.0,
                init_state_concentration=1.0
            )
        )

        # ------------------------------ Fit HMM ------------------------------------ #
        if "DP" in mparams["model_type"]:
            mparams["model_args"].update(K_max=n_states)
        else:
            mparams["model_args"].update(K=n_states)

        # Set seed.
        print "Setting seed to {}.".format(mparams["seed"])
        np.random.seed(mparams["seed"])
        print "Running subject {}, session {}.".format(subj_sess[i_dataset][0], subj_sess[i_dataset][1])

        # Set duration distribution/parameters for semi-Markovian models.
        if "HSMM" in mparams["model_type"]:
            if "Duration" in mparams["model_type"]:
                # Parameters will be used to create duration distribution in constructor
                # method for classes explicitly named after a distribution, e.g.
                # PoissonHSMMPoissonDuration.
                mparams["model_args"].update(alpha_dur=10, beta_dur=2)
            else:
                # Else, if using the general PoissonHSMM class, we will make the
                # duration distributions here first.
                mparams.update(alpha_dur=30, beta_dur=2)
                mparams["model_args"].update(dur_distns \
                                                 =[PoissonDuration(alpha_0=mparams["alpha_dur"], beta_0=mparams["beta_dur"])
                                                   for state in range(n_states)])

        # Instantiate model.
        if mparams["model_type"] == "PoissonHMM":
            hmm = pyhsmm_spiketrains.models.PoissonHMM(**mparams["model_args"])
        elif mparams["model_type"] == "PoissonHSMM":
            hmm = pyhsmm_spiketrains.models.PoissonHSMM(**mparams["model_args"])
        elif mparams["model_type"] == "PoissonHSMMPoissonDuration":
            hmm = pyhsmm_spiketrains.models.PoissonHSMMPoissonDuration(**mparams["model_args"])
        elif mparams["model_type"] == "PoissonHDPHSMM":
            hmm = pyhsmm_spiketrains.models.PoissonHDPHSMM(**mparams["model_args"])

        # Add data.
        hmm.add_data(neural_data)

        # Fit the model with Gibbs sampling
        for itr in progprint_xrange(mparams["n_iter"]):
            hmm.resample_model()
            hmm.resample_obs_hypers_hmc()

        # ---------------------------- Postprocessing ------------------------------- #
        # Get the inferred state sequence
        used_states = hmm.used_states
        n_used_states = len(list(hmm.used_states))
        print("{} states used.".format(n_used_states))

        # Generate state sequence and get firing rate matrix (tuning curves).
        hidden_states = hmm.stateseqs[0]
        fr_mat = hmm.rates

        # Posteriors.
        posteriors = hmm.heldout_state_marginals(neural_data)

        # Get parameters based on model.
        if "HSMM" in mparams["model_type"]:
            trans_mat = hmm.trans_distn.full_trans_matrix
            trans_mat_used = hmm.trans_distn.trans_matrix
            lambdas = [hmm.dur_distns[i_distn].lmbda for i_distn in range(n_used_states)]
        else:
            trans_mat = hmm.trans_distn.trans_matrix
            lambdas = []
        durations = hmm.durations

        # Apply hungarian algorithm if desired.
        APPLY_HUNGARIAN = False
        if APPLY_HUNGARIAN:
            # Apply algorithm and get assignments and overlap.
            row_ind, col_ind, overlap, cost \
                = get_optimal_assignment(labels, hidden_states, method="cosine_sim")
            occupancy = get_overlap_mat(labels, hidden_states, method="count")
            occupancy_normed = occupancy / np.sum(occupancy, axis=1)[:, np.newaxis]

            # Reorder based on optimal assignments------------------------------------------------------------------------------

            # For a total of m hidden state labels and n ground truth labels, if m > n, extra_ind has its elements the column
            # indices (in e.g. the overlap matrix) of the m - n hidden states not matched to a ground truth label. Should be
            # empty if m <= n, since every hidden state will be matched to a ground truth state.
            extra_ind = np.setdiff1d(np.arange(overlap.shape[1]), col_ind)

            # Concatenate col_ind and extra_ind to get col_ind_full, which has as its first n elements the n column indices (
            # of the overlap matrix) of the n hidden states identified as the best matches to the ground truth labels (
            # col_ind). The remaining m - n elements are simply the remaining column indices of the overlap matrix (
            # extra_ind). If m <= n, col_ind_full should be equivalent to col_ind, since extra_ind is empty.
            col_ind_full = np.concatenate((col_ind, extra_ind))

            # Reorder columns of matrices.
            overlap_reordered = overlap[:, col_ind_full]
            cost_reordered = cost[:, col_ind_full]
            occupancy_reordered = occupancy[:, col_ind_full]
            occupancy_normed_reordered = occupancy_normed[:, col_ind_full]

            # Reorder transition matrix.
            trans_mat_reordered = trans_mat[:, col_ind_full]
            trans_mat_reordered = trans_mat_reordered[col_ind_full, :]
            trans_mat_used_reordered = trans_mat_used[:, col_ind_full]
            trans_mat_used_reordered = trans_mat_used[col_ind_full, :]

            # Reorder hidden states.
            hidden_states_reordered = np.full(len(hidden_states), np.nan)
            i_extra = 0
            for i_state, state in enumerate(np.unique(hidden_states)):
                # If current hidden state was one of the ones matched to a ground truth state (will be true if m <= n, for m and
                # n as defined above; will be false for some hidden states if m > n), replace that hidden state's label with the
                # ground truth label. Else, retain the original hidden state label.
                if np.isin(state, col_ind):
                    hidden_states_reordered[hidden_states == state] = np.where(col_ind == state)[0][0]
                else:
                    # For hidden states not matched, use labels left over (i.e., not used in col_ind).
                    hidden_states_reordered[hidden_states == state] = extra_ind[i_extra]
                    i_extra += 1
        else:
            row_ind = []
            col_ind = []
            overlap = []
            occupancy = []
            occupancy_normed = []
            cost = []
            extra_ind = []
            col_ind_full = []
            overlap_reordered = []
            cost_reordered = []
            occupancy_reordered = []
            occupancy_normed_reordered = []
            trans_mat_reordered = []
            trans_mat_used_reordered = []
            hidden_states_reordered = []

        # %%
        # ---------------------------------- Save ----------------------------------- #
        if SAVE:
            # Change to save directory and prepare filename.
            os.chdir(save_dir)
            save_filename = load_filename + "-" + mparams["model_type"] + "-" \
                            + "gibbs{}-{}".format(mparams["n_iter"], run_id)
            if data_type == "sim":
                save_filename = save_filename + "-sim"

            # Save postprocessing vars for matlab.
            scipy.io.savemat(save_filename + ".mat",
                             {"used_states": used_states,
                              "hidden_states": hidden_states,
                              "hidden_states_reordered": hidden_states_reordered,
                              "overlap": overlap,
                              "overlap_reordered": overlap_reordered,
                              "occupancy": occupancy,
                              "occupancy_reordered": occupancy_reordered,
                              "cost": cost,
                              "cost_reordered": cost_reordered,
                              "n_used_states": n_used_states,
                              "lambdas": lambdas,
                              "durations": durations,
                              "fr_mat": fr_mat,
                              "trans_mat": trans_mat,
                              "trans_mat_used": trans_mat_used,
                              "trans_mat_reordered": trans_mat_reordered,
                              "trans_mat_used_reordered": trans_mat_used_reordered,
                              "col_ind": col_ind,
                              "row_ind": row_ind,
                              "extra_ind": extra_ind,
                              "col_ind_full": col_ind_full,
                              "posteriors": posteriors,
                              "pparams": pparams,
                              "mparams": mparams})

            # Save postprocessing vars for Python.
            with open(save_filename, "wb") as file_:
                pickle.dump([hmm,
                             pparams,
                             mparams,
                             row_ind,
                             col_ind,
                             extra_ind,
                             col_ind_full,
                             used_states,
                             hidden_states,
                             hidden_states_reordered,
                             overlap,
                             cost,
                             occupancy,
                             overlap_reordered,
                             cost_reordered,
                             occupancy_reordered,
                             n_used_states,
                             lambdas,
                             durations,
                             fr_mat,
                             trans_mat,
                             trans_mat_used,
                             trans_mat_reordered,
                             trans_mat_used_reordered,
                             posteriors],
                            file_)

    os.chdir(curr_dir)
