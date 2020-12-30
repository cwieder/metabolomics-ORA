# Tools for pathway analysis in metabolomics

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import molmass
import time

pd.options.mode.chained_assignment = None  # default='warn'

def data_processing(raw_matrix, firstrow, firstcol):
    '''
    Filtering low abundance metabolites, data cleaning and imputation using minimum value.
    :param raw_matrix: raw abundance matrix with n-samples and m-metabolites
    :param firstrow: First row containing values (integer)
    :return: imputed, log-transformed and standardised matrix
    '''

    # Remove commas and convert to numeric
    processed_matrix = raw_matrix.iloc[firstrow:, firstcol:]
    processed_matrix = processed_matrix.replace(',', '', regex=True)
    processed_matrix = processed_matrix.apply(pd.to_numeric)
    processed_matrix.replace(0, np.nan, inplace=True)

    processed_matrix = processed_matrix.loc[:, processed_matrix.isnull().mean() < 0.9]
    # Remove metabolites not present in > 90% of samples
    # by indexing df by rows with metabolites present in more than 90% of samples

    # Missing value imputation using minimum value/2
    imputed_matrix = processed_matrix.replace(np.nan, processed_matrix.min(axis=0)/2)
    # Log2 transformation
    log2_matrix = np.log2(imputed_matrix)
    # Standardisation by mean centering and scaling to unit variance
    standarised_mat = StandardScaler().fit_transform(log2_matrix)
    log2_matrix.loc[:, :] = standarised_mat
    return log2_matrix

def plot_PCA(matrix, metadata, title, n_comp=100):
    '''
    PCA plot
    :param matrix: processed data matrix
    :param metadata: pandas series with group names
    :param title: plot title
    :param n_comp: number of PCA components
    :return: PCA plot
    '''

    pca = PCA(n_components=n_comp)
    projected = pca.fit_transform(matrix)
    print(round(sum(list(pca.explained_variance_ratio_))*100, 2), "%")

    scatter_x = projected[:,0]
    scatter_y = projected[:,1]
    samples = matrix.index.tolist()
    group = metadata
    uniq_sample = np.unique(group)
    cmap = ['tab:green', 'tab:orange', 'tab:blue', 'red', 'cyan', 'magenta', 'brown', 'yellow']
    cdict = {samp: cmap[num] for num, samp in enumerate(uniq_sample)}

    plt.style.use("ggplot")
    fig, ax = plt.subplots()
    for g in np.unique(group):
        ix = np.where(group == g)
        ax.scatter(scatter_x[ix], scatter_y[ix], c=cdict[g], label=g, s=15)
    ax.legend()

    plt.title("PCA for " + title)
    plt.xlabel("Component 1:  " + str(round(pca.explained_variance_ratio_[0]*100, 2)) + "%")
    plt.ylabel("Component 2: " + str(round(pca.explained_variance_ratio_[1]*100, 2)) + "%")
    plt.savefig(title + ".png")
    plt.show()

def linear_regression(matrix, metadatadict):
    # Add new columns based on metadata from dict
    matrix['PMH_status'] = matrix.index.map(lambda x: metadatadict[x][0])
    matrix['Target'] = pd.factorize(matrix['PMH_status'])[0]
    coefs = []
    pvals = []
    metabolites = matrix.columns.tolist()[:-2]
    for metabolite in metabolites:
        X = matrix[str(metabolite)]
        y = matrix['Target']
        model = sm.OLS(y, X.astype(float)).fit()
        coefs.append(model.params[0])
        pvals.append(model.pvalues[0])
    padj = sm.stats.multipletests(pvals, 0.05, method="bonferroni")
    print("Corrected alpha:", padj[3])
    results = pd.DataFrame(zip(metabolites, coefs, pvals, padj[1]), columns=["Metabolite", "Fold-change", "P-value", "P-adjust"])
    return results

def t_tests(matrix, classes, multiple_correction_method):

    matrix['Target'] = pd.factorize(classes)[0]
    metabolites = matrix.columns.tolist()[:-1]

    pvalues = []
    for metabolite in metabolites:
        # print(matrix[matrix['Target'] == 0][metabolite].tolist())
        if isinstance(matrix[matrix['Target'] == 0][metabolite], pd.DataFrame):
            print(metabolite)
            print("duplicate compounds found")
            print(matrix[matrix['Target'] == 0][metabolite])
        group1 = matrix[matrix['Target'] == 0][metabolite].tolist()
        group2 = matrix[matrix['Target'] == 1][metabolite].tolist()
        stat, pval = stats.ttest_ind(group1, group2)
        pvalues.append(pval)
    padj = sm.stats.multipletests(pvalues, 0.05, method=multiple_correction_method)
    results = pd.DataFrame(zip(metabolites, pvalues, padj[1]),
                           columns=["Metabolite", "P-value", "P-adjust"])
    return results

def over_representation_analysis(DEM_list, background_list, pathways_df):
    """
    Function for over representation analysis
    :param DEM_list: List of differentially exprssed metabolite IDENTIFIERS
    :param background_list: background list of IDENTIFIERS
    :param pathways_df: pathway dataframe containing compound identifiers
    :return: DataFrame of ORA results for each pathway, p-value, q-value, hits ratio
    """
    KEGG_pathways = pathways_df.dropna(axis=0, how='all', subset=pathways_df.columns.tolist()[1:])
    pathways = KEGG_pathways.index.tolist()
    pathway_names = KEGG_pathways["Pathway_name"].tolist()
    pathway_dict = dict(zip(pathways, pathway_names))
    KEGG_pathways.drop('Pathway_name', axis=1, inplace=True)

    pathways_with_compounds = []
    pathway_names_with_compounds = []
    pvalues = []
    pathway_ratio = []
    for pathway in pathways:
        # perform ORA for each pathway
        pathway_compounds = KEGG_pathways.loc[pathway, :].tolist()
        pathway_compounds = [i for i in pathway_compounds if str(i) != "nan"]
        if not pathway_compounds or len(pathway_compounds) < 3:
            # ignore pathway if contains no compounds or has less than 3 compounds
            continue
        else:
            DEM_in_pathway = len(set(DEM_list) & set(pathway_compounds))
            # k: compounds in DEM list AND pathway
            DEM_not_in_pathway = len(np.setdiff1d(DEM_list, pathway_compounds))
            # K: compounds in DEM list not in pathway
            compound_in_pathway_not_DEM = len(set(np.setdiff1d(background_list, DEM_list)) & set(pathway_compounds))
            # compounds present in bg set and pathway
            compound_not_in_pathway_not_DEM = len(np.setdiff1d(np.setdiff1d(background_list, DEM_list), pathway_compounds))
            # compounds in background list not present in pathway

            if (DEM_in_pathway and compound_in_pathway_not_DEM) == 0:
                # ignore pathway if there are no DEM/background list compounds in that pathway
                continue
            else:
                # Create 2 by 2 contingency table
                pathway_ratio.append(str(DEM_in_pathway) + "/" + str(len(pathway_compounds)))
                pathways_with_compounds.append(pathway)
                pathway_names_with_compounds.append(pathway_dict[pathway])
                contingency_table = np.array([[DEM_in_pathway, compound_in_pathway_not_DEM],
                                              [DEM_not_in_pathway, compound_not_in_pathway_not_DEM]])

                # Run right tailed Fisher's exact test
                oddsratio, pvalue = stats.fisher_exact(contingency_table, alternative="greater")
                pvalues.append(pvalue)
    try:
        padj = sm.stats.multipletests(pvalues, 0.05, method="fdr_bh")
        results = pd.DataFrame(zip(pathways_with_compounds, pathway_names_with_compounds, pathway_ratio, pvalues, padj[1]),
                               columns=["Pathway_ID", "Pathway_name", "Hits", "P-value", "P-adjust"])
    except ZeroDivisionError:
        padj = [1] * len(pvalues)
        results = pd.DataFrame(zip(pathways_with_compounds, pathway_names_with_compounds, pathway_ratio, pvalues, padj),
                           columns=["Pathway_ID", "Pathway_name", "Hits", "P-value", "P-adjust"])
    return results


def reduce_background_list_ora(background_list, percentage, DEM_list, pathways_df, keep_DEM=False):
    """
    Reduces size of background list by random removal of compounds
    :param background_list: background list of compound names/IDs
    :param percentage: percentage reduction of list desired
    :param pathways_df: df of organism specific pathways
    :param keep_DEM: do not delete differentially abundant metabolites
    :return: reduced background list
    """
    bg_list_without_dem = [i for i in background_list if i not in DEM_list]
    if not keep_DEM:
        list_size = int(len(background_list)*(percentage/100))
    else:
        list_size = int(len(bg_list_without_dem)*(percentage/100))
    print(list_size, str(percentage) +"%")

    p_vals = []
    q_vals = []
    baseline_significant_paths = []
    proportion_of_original_pathways_signficant_p = []

    # proportion_of_original_pathways_signficant_q = []

    baseline_res = over_representation_analysis(DEM_list, background_list, pathways_df)
    baseline_significant_paths.append(len(baseline_res[baseline_res["P-value"] < 0.1]["P-value"].tolist()))
    # q_vals.append(len(baseline_res[baseline_res["P-adjust"] < 0.1]["P-adjust"].tolist()))

    for i in range(0, 100):
        if not keep_DEM:
            bg_list_reduced = np.random.choice(background_list, list_size, replace=False)
        else:
            bg_list_reduced = np.random.choice(np.setdiff1d(bg_list_without_dem, DEM_list), list_size, replace=False)
            bg_list_reduced = bg_list_reduced.tolist() + DEM_list
        ora_res = over_representation_analysis(DEM_list, bg_list_reduced, pathways_df)
        p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
        print(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()), baseline_significant_paths[0])
        proportion_of_original_pathways_signficant_p.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist())/baseline_significant_paths[0])
        q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
        # proportion_of_original_pathways_signficant_p.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist())/q_vals[0])
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    mean_proportion_p_vals = np.mean(proportion_of_original_pathways_signficant_p)
    sd_p_signficant_paths = stats.sem(p_vals)
    sd_q_signficant_paths = stats.sem(q_vals)
    sd_proportion_p_vals = stats.sem(proportion_of_original_pathways_signficant_p)
    return [mean_p_signficant_paths, mean_q_signficant_paths, mean_proportion_p_vals,
            sd_p_signficant_paths, sd_q_signficant_paths, sd_proportion_p_vals]

def misidentify_metabolites(percentage, processed_matrix, organism_compounds, background_list, pathway_df,
                            zamboni=False):
    '''
    Randomly swaps a percentage of KEGG compounds and then performs ORA
    :param percentage: percentage of compounds to be misidentified
    :param processed_matrix: processed abundance matrix with KEGG compounds as columns
    :param organism_compounds: all KEGG compounds for the organism
    :param background_list: background list for dataset
    :param pathway_df: list of KEGG pathways
    :param zamboni: False, assign true if dataset is Zamboni
    :return: mean number of p-values significant at P <0.1, Q-values, and standard deviation
    '''

    p_vals = []
    q_vals = []
    significant_pathways = [] # at P < 0.1

    if zamboni == False:
        mat_unannotated = processed_matrix.iloc[:, :-1]
        metabolites = mat_unannotated.columns.tolist()
        n_misidentified = int(len(metabolites) * (percentage / 100))
        for i in range(0, 100):
            # Randomly replace n compounds
            metabolites_to_replace = np.random.choice(metabolites, n_misidentified, replace=False)
            replacement_compounds = np.random.choice(np.setdiff1d(organism_compounds, background_list), n_misidentified, replace=False)
            # TODO: CHECK THAT ORGANISM COMPOUNDS CONTAIN THE RIGHT THINGS e.g. glycans? drugs? should be included?
            replacement_dict = dict(zip(metabolites_to_replace, replacement_compounds))
            misidentified_matrix = mat_unannotated.rename(columns=replacement_dict)

            # Perform t-tests and ORA
            ttest_res = t_tests(misidentified_matrix, processed_matrix["Group"], "fdr_bh")
            DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
            significant_pathways.append(ora_res[ora_res["P-value"] < 0.1]["Pathway_ID"].tolist())

    elif zamboni == True:
        metabolites = processed_matrix.columns.tolist()
        n_misidentified = int(len(set(metabolites)) * (percentage / 100))
        for i in range(0, 100):
            metabolites_to_replace = np.random.choice(list(set(metabolites)), n_misidentified, replace=False)
            replacement_compounds = np.random.choice(np.setdiff1d(organism_compounds,
                                                                  [background_list + list(metabolites_to_replace)]), n_misidentified,
                                                     replace=False)
            replacement_dict = dict(zip(metabolites_to_replace, replacement_compounds))
            misidentified_matrix = processed_matrix.rename(columns=replacement_dict)
            DEM = []
            for x in misidentified_matrix.T.itertuples():
                if x[1] > 6 or x[1] < -6:
                    DEM.append(x[0])
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
            significant_pathways.append(ora_res[ora_res["P-value"] < 0.1]["Pathway_ID"].tolist())
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    sd_p_signficant_paths = stats.sem(p_vals)
    sd_q_signficant_paths = stats.sem(q_vals)

    return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths, significant_pathways]

def misidentify_metabolites_by_mass(percentage, processed_matrix, pathway_df, all_compound_masses, organism_bg, zamboni=False):
    '''
    Randomly swaps a percentage of KEGG compounds and then performs ORA
    :param percentage: percentage of compounds to be misidentified
    :param processed_matrix: processed abundance matrix with KEGG compounds as columns
    :param pathway_df: list of KEGG pathways
    :param all_compound_masses: list of all compounds and masses
    :param organism_bg: organism specific compounds
    :return: mean number of p-values significant at P <0.1, Q-values, and standard deviation
    '''

    p_vals = []
    q_vals = []
    KEGG_compounds_masses_organism = all_compound_masses[all_compound_masses.index.isin(organism_bg)]

    if not zamboni:
        mat_unannotated = processed_matrix.iloc[:, :-1]
        metabolites = list(set(mat_unannotated.columns.tolist()))
        print(len(metabolites))
        n_misidentified = int(len(set(metabolites)) * (percentage / 100))

        mass_dict = dict(zip(KEGG_compounds_masses_organism.index.tolist(), KEGG_compounds_masses_organism['molecular_weight'].tolist()))
        metabolites_mass_dict = {k: v for k, v in mass_dict.items() if k in metabolites}
        misidentifiable_metabolites = dict()

        for cpd, mass in metabolites_mass_dict.items():
            # mass_window = (mass - 5, mass + 5)
            mass_window = ((mass - (mass/1000000)*20), (mass + (mass/1000000)*20))
            cpd_info = KEGG_compounds_masses_organism[KEGG_compounds_masses_organism['molecular_weight'].between(mass_window[0], mass_window[1],
                                                                inclusive=False)].index.tolist()
            if len(cpd_info) > 1:
                misidentifiable_metabolites[cpd] = np.setdiff1d(cpd_info, cpd).tolist()
        print(len(misidentifiable_metabolites))
        for i in range(0, 100):
            replacement_dict = dict()
            while len(replacement_dict) < n_misidentified:
                cpd_to_relpace = np.random.choice(list(misidentifiable_metabolites.keys()), 1)[0]
                replacement_cpd = np.random.choice(misidentifiable_metabolites[cpd_to_relpace], 1)[0]
                if replacement_cpd not in list(replacement_dict.values()) and replacement_cpd not in metabolites:
                    replacement_dict[cpd_to_relpace] = replacement_cpd
            print(replacement_dict)
            misidentified_matrix = mat_unannotated.rename(columns=replacement_dict)
            # Perform t-tests and ORA
            ttest_res = t_tests(misidentified_matrix, processed_matrix["Group"], "fdr_bh")
            DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    elif zamboni:
        metabolites = list(set(processed_matrix.columns.tolist()))
        print(len(metabolites))
        n_misidentified = int(len(set(metabolites)) * (percentage / 100))

        mass_dict = dict(zip(KEGG_compounds_masses_organism.index.tolist(), KEGG_compounds_masses_organism['molecular_weight'].tolist()))
        metabolites_mass_dict = {k: v for k, v in mass_dict.items() if k in metabolites}
        misidentifiable_metabolites = dict()
        for cpd, mass in metabolites_mass_dict.items():
            # mass_window = (mass - 0.004, mass + 0.004)
            mass_window = ((mass - (mass / 1000000) * 20), (mass + (mass / 1000000) * 20))
            cpd_info = KEGG_compounds_masses_organism[KEGG_compounds_masses_organism['molecular_weight'].between(mass_window[0], mass_window[1],
                                                                inclusive=False)].index.tolist()
            if len(cpd_info) > 1:
                misidentifiable_metabolites[cpd] = np.setdiff1d(cpd_info, cpd).tolist()

        print(len(misidentifiable_metabolites))
        for i in range(0, 100):
            replacement_dict = dict()
            while len(replacement_dict) < n_misidentified:
                cpd_to_relpace = np.random.choice(list(misidentifiable_metabolites.keys()), 1)[0]
                replacement_cpd = np.random.choice(misidentifiable_metabolites[cpd_to_relpace], 1)[0]
                if replacement_cpd not in list(replacement_dict.values()) and replacement_cpd not in metabolites and cpd_to_relpace not in list(replacement_dict.keys()):
                    replacement_dict[cpd_to_relpace] = replacement_cpd
            misidentified_matrix = processed_matrix.rename(columns=replacement_dict)
            DEM = []
            for x in misidentified_matrix.T.itertuples():
                if x[1] > 6 or x[1] < -6:
                    DEM.append(x[0])
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    sd_p_signficant_paths = stats.sem(p_vals)
    sd_q_signficant_paths = stats.sem(q_vals)
    return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths]

def misidentify_metabolites_by_formula(percentage, processed_matrix, pathway_df, all_cpd_formulas, organism_bg,
                                       zamboni=False):
    '''
    Randomly swaps a percentage of KEGG compounds and then performs ORA
    :param percentage: percentage of compounds to be misidentified
    :param processed_matrix: processed abundance matrix with KEGG compounds as columns
    :param pathway_df: list of KEGG pathways
    :return: mean number of p-values significant at P <0.1, Q-values, and standard deviation
    '''
    KEGG_compounds_formula_organism = all_cpd_formulas[all_cpd_formulas.index.isin(organism_bg)]
    print(len(organism_bg))
    print(len(KEGG_compounds_formula_organism))
    p_vals = []
    q_vals = []

    if not zamboni:
        mat_unannotated = processed_matrix.iloc[:, :-1]
        metabolites = mat_unannotated.columns.tolist()
        n_misidentified = int(len(metabolites)*(percentage/100))

        misidentifiable_metabolites = dict()
        compound_formula_dict = dict(zip(KEGG_compounds_formula_organism.index.tolist(), KEGG_compounds_formula_organism['molecular_formula'].tolist()))
        metabolites_formula_dict = {k: v for k, v in compound_formula_dict.items() if k in metabolites}
        # if at least 1 other formula, add to misidentifiable metabolites dict

        for cpd, formula in metabolites_formula_dict.items():
            cpd_info = KEGG_compounds_formula_organism[
                KEGG_compounds_formula_organism['molecular_formula'] == formula].index.tolist()
            if len(cpd_info) > 1:
                misidentifiable_metabolites[cpd] = np.setdiff1d(cpd_info, cpd).tolist()
        print(len(misidentifiable_metabolites))
        for i in range(0, 100):
            replacement_dict = dict()
            while len(replacement_dict) < n_misidentified:
                cpd_to_relpace = np.random.choice(list(misidentifiable_metabolites.keys()), 1)[0]
                replacement_cpd = np.random.choice(misidentifiable_metabolites[cpd_to_relpace], 1)[0]
                if replacement_cpd not in list(replacement_dict.values()) and replacement_cpd not in metabolites:
                    replacement_dict[cpd_to_relpace] = replacement_cpd
            print(replacement_dict)
            misidentified_matrix = mat_unannotated.rename(columns=replacement_dict)
            # Perform t-tests and ORA
            ttest_res = t_tests(misidentified_matrix, processed_matrix["Group"], "fdr_bh")
            DEM = ttest_res[ttest_res["P-adjust"] < 0.05]["Metabolite"].tolist()
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    elif zamboni:
        metabolites = list(set(processed_matrix.columns.tolist()))
        print(len(metabolites))
        n_misidentified = int(len(set(metabolites)) * (percentage / 100))
        misidentifiable_metabolites = dict()
        compound_formula_dict = dict(zip(KEGG_compounds_formula_organism.index.tolist(),
                                         KEGG_compounds_formula_organism['molecular_formula'].tolist()))
        metabolites_formula_dict = {k: v for k, v in compound_formula_dict.items() if k in metabolites}
        # if at least 1 other formula, add to misidentifiable metabolites dict

        for cpd, formula in metabolites_formula_dict.items():
            cpd_info = KEGG_compounds_formula_organism[
                KEGG_compounds_formula_organism['molecular_formula'] == formula].index.tolist()
            if len(cpd_info) > 1:
                misidentifiable_metabolites[cpd] = np.setdiff1d(cpd_info, cpd).tolist()
        print(len(misidentifiable_metabolites))
        for i in range(0, 100):
            replacement_dict = dict()
            while len(replacement_dict) < n_misidentified:
                cpd_to_relpace = np.random.choice(list(misidentifiable_metabolites.keys()), 1)[0]
                replacement_cpd = np.random.choice(misidentifiable_metabolites[cpd_to_relpace], 1)[0]
                if replacement_cpd not in list(replacement_dict.values()) and replacement_cpd not in metabolites:
                    replacement_dict[cpd_to_relpace] = replacement_cpd
            print(replacement_dict)
            misidentified_matrix = processed_matrix.rename(columns=replacement_dict)
            DEM = []
            for x in misidentified_matrix.T.itertuples():
                if x[1] > 6 or x[1] < -6:
                    DEM.append(x[0])
            ora_res = over_representation_analysis(DEM, misidentified_matrix.columns.tolist(), pathway_df)
            p_vals.append(len(ora_res[ora_res["P-value"] < 0.1]["P-value"].tolist()))
            q_vals.append(len(ora_res[ora_res["P-adjust"] < 0.1]["P-adjust"].tolist()))
    mean_p_signficant_paths = np.mean(p_vals)
    mean_q_signficant_paths = np.mean(q_vals)
    sd_p_signficant_paths = stats.sem(p_vals)
    sd_q_signficant_paths = stats.sem(q_vals)
    return [mean_p_signficant_paths, mean_q_signficant_paths, sd_p_signficant_paths, sd_q_signficant_paths]