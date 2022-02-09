"""
Direction prediction based on learning dataset from reactome
PPI direction calculated from domain interaction directions
"""

# Imports
import sqlite3, csv, os
import pandas as pd
import logging
import pickle
# # Initiating logger
# logger = logging.getLogger()
# handler = logging.FileHandler('../../workflow/SLK3.log')
# logger.setLevel(logging.DEBUG)
# handler.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# handler.setFormatter(formatter)
# logger.addHandler(handler)


class DirScore:
    def __init__(self):
        # Defining constants
        self.REACTOME_DB = '../../ARNlib/mapper/protein/output/reactome_mapped.db'
        self.PFAM_FILE = ['../prediction/direction/files/uniprot-pfam_human.tab']
        logging.basicConfig(level=logging.DEBUG)
        self.pfam_dict = {}
        self.dir_score_dict = {}
        # Adding the two output dictionaries of test_score function to a pickle files
        # so that the next function can access them inbetween script executions
        # TODO: remove pickle files after each run
        self.PICKLE_FILE = 'dir_score.pickle'
        if os.path.isfile(self.PICKLE_FILE):
            self.pfam_dict, self.dir_score_dict = pickle.load(open(self.PICKLE_FILE, 'rb'))
        else:
            self.test_scores()
            pickle.dump((self.pfam_dict, self.dir_score_dict), open(self.PICKLE_FILE, 'wb'))

    def test_scores(self):
            # Setting as global so next script can access it
            df_all = pd.DataFrame(columns=['a_dom', 'b_dom'])

            conn = sqlite3.connect(self.REACTOME_DB)
            # Setting up learning data set
            logging.debug("Started connection to reactome dataset")
            for inpfam in self.PFAM_FILE:
                with open(inpfam) as infile:
                    infile.readline()
                    for line in infile:
                        line = line.strip().split('\t')
                        if len(line) == 4:
                            self.pfam_dict[line[0]] = line[3].split(';')[0:-1]

            with conn:

                c = conn.cursor()

                counter = 0
                # Getting PPI data
                logging.debug('Getting PPI data')
                c.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM edge")
                while True:
                    row = c.fetchone()
                    counter += 1

                    if row is None:
                        break
                    else:
                        a_node = row[0].split(':')[1]
                        b_node = row[1].split(':')[1]
                        if a_node not in self.pfam_dict or b_node not in self.pfam_dict:
                            continue

                        int_list = [self.pfam_dict[a_node], self.pfam_dict[b_node]]

                        for id1, id2 in zip(int_list[0], int_list[1]):
                            # Setting up dataframe for all domain-domain interactions
                            # len(df_all) sets the name of the line
                            # df_all = df_all.loc(len(df_all), col=['a_dom', 'b_dom'], value=[id1, id2])
                            df_all.loc[len(df_all), 'a_dom'] = id1
                            df_all.loc[len(df_all), 'b_dom'] = id2

                            # All domains in a dataframe, without direction
            all_domain_df = df_all['a_dom'].append(df_all['b_dom']).reset_index(name='domain')
            all_count = all_domain_df.groupby('domain').size().reset_index(name='counter')

            # Getting probability of each domain
            # Number of domain occurrence / Number of all domains
            logging.debug('Getting probability of each domain')
            prob_dom = {}
            # Number of all domain occurrences
            total_occurrence = all_count['counter'].sum()
            # Iterating over domains
            for index, domain in all_count['domain'].iteritems():
                dom_count = all_count.loc[all_count['domain'] == domain, 'counter'].iloc[0]
                P_domain = dom_count / total_occurrence
                # Adding data into a dictionary
                prob_dom[domain] = P_domain
                #print(domain, P_domain)

            # Getting directed domain-domain interaction probabilities
            # Number of directed DDI / number of all DDIs
            logging.debug('Getting DDI probabilities')
            prob_inter = {}
            # Getting the occurrences for each directed interaction
            all_inter_counted = df_all.groupby(['a_dom', 'b_dom']).size().reset_index(name='counter')
            all_inter_counter = all_inter_counted['counter'].sum()
            # Iterating over interactions
            for index2, count in all_inter_counted['counter'].iteritems():
                P_inter = count / all_inter_counter
                # Getting domain ids
                a_dom = all_inter_counted.loc[all_inter_counted['counter'] == count, 'a_dom'].iloc[0]
                b_dom = all_inter_counted.loc[all_inter_counted['counter'] == count, 'b_dom'].iloc[0]
                # Adding the into a dictionary
                prob_inter['->'.join((a_dom, b_dom))] = P_inter

            # Calculating direction score
            # (P_AtoB - P_BtoA) / P_A * P_B
            logging.debug('Calculating direction scores')

            for key in prob_inter.keys():
                a = key.split('->')[0]
                b = key.split('->')[1]
                other_dir = '->'.join((b, a))
                if other_dir in prob_inter.keys():
                    dir_score = (prob_inter[key] - prob_inter[other_dir]) / prob_dom[a] * prob_dom[b]
                    self.dir_score_dict[key] = dir_score
                else:
                    dir_score = (prob_inter[key] - 0) / prob_dom[a] * prob_dom[b]
                    self.dir_score_dict[key] = dir_score
                    #print(key, dir_score)

            #return self.dir_score_dict, self.pfam_dict


    # LAYER 3
    def apply_to_db(self):
        #logger.debug(self.pfam_dict)
        #logger.debug(self.dir_score_dict)
        conn2 = sqlite3.connect('ARN2_layers.db')
        # logger.debug("Connected to '%s" % conn2)
        with conn2:
            c2 = conn2.cursor()
            c22 = conn2.cursor()
            c2.execute("SELECT interactor_a_node_name, interactor_b_node_name FROM layer3")
            while True:
                row = c2.fetchone()
                if row is None:
                    break
                else:
                    prot_a = row[0].split(':')[1]
                    prot_b = row[1].split(':')[1]
                    dir_score_sum = 0
                    # Summing DDI scores
                    #logging.debug('Summing DDI scores')
                    if prot_a in self.pfam_dict.keys() and prot_b in self.pfam_dict.keys():
                        for dom_a, dom_b in zip(self.pfam_dict[prot_a], self.pfam_dict[prot_b]):
                            #print(dir_score_dict['->'.join((dom_a, dom_b))])
                            if '->'.join((dom_a, dom_b)) in self.dir_score_dict.keys():
                                dir_score_sum += self.dir_score_dict['->'.join((dom_a, dom_b))]
                        # To get final direction score of the unknown PPIs we calculate
                        # the average of each proteins' all domain interaction scores
                        if len(self.pfam_dict[prot_a]) * len(self.pfam_dict[prot_b]) == 0:
                            logging.debug(prot_a, len(self.pfam_dict[prot_a]), prot_b, len(self.pfam_dict[prot_b]))
                            continue
                        else:
                            dir_score_final_PPI = dir_score_sum / (len(self.pfam_dict[prot_a]) * len(self.pfam_dict[prot_b]))
                        #logging.debug("Updating scores")
                        c22.execute("UPDATE layer3 SET confidence_scores = '%s' "
                                    "WHERE layer3.interactor_a_node_name = '%s' AND layer3.interactor_b_node_name = '%s'"
                                    % ('|dir_pred:' + str(dir_score_final_PPI), row[0], row[1]))


if __name__ == '__main__':
    test = DirScore()
    logger.debug('Creating test set')
    test.test_scores()
    logger.debug('Adding scores to dataset')
    test.apply_to_db()
    logger.debug('Direction prediction done')





