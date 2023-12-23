import sys
import random
import copy

sequence_file = sys.argv[1]
bootstraps = 100

# open file from command line argument
f = open(sequence_file, "r");

sequence_list = []
sequence_ids = []

# getting the sequences / their ids into two lists
# sequence and the id will have the same index in their respective list
while True:
    the_sequence = f.readline()
    the_sequence = the_sequence.rstrip()
    if not the_sequence:
        break
    if(the_sequence.startswith(">")):
        sequence_ids.append(the_sequence[1:])
    else:
        sequence_list.append(the_sequence)

f.close()


# this assumes the lengths of sequences are the same
# finds number of differences in sequences
# then divides by length to get p distance
def compareSeqs(seq1, seq2, length):
    differences = 0
    for i in range(0, length):
        if(seq1[i] != seq2[i]):
            differences += 1

    return differences/length


# lengths of all sequences are the same
sequence_length = len(sequence_list[0])

# this gives the genetic distance matrix
# having it as a function helps with bootstraps 
def get_distances(sequenceL):

    # square matrix, each id has a corresponding sequence
    distance_matrix = [[0.0 for _ in range(len(sequenceL))] for _ in range(len(sequenceL))]


    # for each item in the sequence list, compare with every other item
    # use compareSeqs to get number of differences/length and round to 15 decimal places like in genetic_distance.txt
    # do this for every pair of sequences except duplicates, so no comparing taxon1 to itself

    for i in range(0, len(sequence_list)):
        for j in range(0, len(sequence_list)):
            if (i != j):
                distance_matrix[i][j] = (compareSeqs(sequenceL[i], sequenceL[j], sequence_length))
    
    return distance_matrix


# stores genetic distance matrix for input file
true_distance_matrix = get_distances(sequence_list)


# write genetic distances to file with tab delimited format

distances_file = open("genetic-distances.txt", "w")

for i in range(0, len(sequence_ids)):
    if(i == len(sequence_ids) - 1):
        distances_file.write(sequence_ids[i])
        continue
    else:
        distances_file.write(sequence_ids[i] + '\t')


for i in range(0, len(sequence_list)):
    distances_file.write('\n')
    for j in range(0, len(sequence_list)):
        if(j == 0):
            distances_file.write(sequence_ids[i] + '\t' + str(true_distance_matrix[i][0]) + '\t')
        elif(j == len(sequence_list) - 1):
            distances_file.write(str(true_distance_matrix[i][j]))
        else:
            distances_file.write(str(true_distance_matrix[i][j]) + '\t')
    

distances_file.close()


# sums designated row of matrix
def sumRow(row, d_matrix):
    sum = 0
    for i in range(0, len(d_matrix[row])):
        sum += d_matrix[row][i]
    return sum

# calculates the Q matrix 
# skips over indices that are equal, like distance_matrix[0][0]
# same size as distance matrix just changed values 
# uses formula (n-2)*d(i, j) - (sum of row i) - (sum of row j)

def calcQMatrix(d_matrix):
    n = len(d_matrix)
    q_matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    for i in range(0, n):
        for j in range(0, n):
            if(i != j):
                q_matrix[i][j] = (n - 2) * (d_matrix[i][j]) - (sumRow(i, d_matrix)) - (sumRow(j, d_matrix))
    
    return q_matrix

# does the distance of each item to the Node
# so if (a, b) makes node c, finds distance from a to c and b to c
# returns two element list with the two distances 
# [distance a to c, distance b to c]
# uses formula (1/2) * d(a, b) + (1/(2*(n-2))) * (sum of row a - sum of row b)

def distanceToNode(index, d_matrix):
    a = index[0]
    b = index[1]
    n = len(d_matrix)
    distances = [0, 0]
    distances[0] = (0.5 * d_matrix[a][b]) + (1/(2*(n-2)))*(sumRow(a, d_matrix) - sumRow(b, d_matrix))
    distances[1] = d_matrix[a][b] - distances[0]

    return distances



# returns first instance of smallest item in matrix
# so if there are two places with the same number as the smallest
# returns the first one
def findMin(d_matrix):
    n = len(d_matrix)
    minIndex = [0, 0]
    # some random number that is larger than any distance
    min = 132423
    for i in range(0, n):
        for j in range(0, n):
            if(d_matrix[i][j] < min and d_matrix[i][j] != min):
                min = d_matrix[i][j]
                minIndex = [i, j]
    
    
    return minIndex


# updates distance matrix after joined nodes

def updateDistanceMatrix(index, d_matrix):
    n = len(d_matrix)
    replace = index[0]
    delete = index[1]
    new_matrix = []

    # updated values are in the rows / columns of the joined node
    # so if the joined nodes were in row/columns [3, 4]
    # 3 will be the row of the new node, 4 will be deleted from matrix
    # only row/column 3 will be changed, other values will be the same as previous distance matrix


    
    # this is just transferring over the previous matrix items to the new matrix
    # basically just moves things to corresponding index, unless the skipped row
    # or column has been covered, then it will just add to index - 1
    # for example if the skipped row has been went over but not the skipped column
    # it will put distance_matrix[row][col] in new_matrix[row - 1][col] 

    for i in range(0, n):
        if(i != delete):
            new_row = []
            for j in range(0, n):
                if(j != delete):
                    if(j == replace or i == replace):
                        new_row.append(0)
                    else:
                        new_row.append(d_matrix[i][j])

            new_matrix.append(new_row)
    
    

    # this is the calculations for new distances in the row/columns
    # if c is the entry and the distance to new node made up of (a, b) is calculated
    # it uses formula d(c, new node) = 0.5 * (d(a, c) + d(b, c) - d(a, b))

    behind_row = 0
    for i in range(0, n):
        if(i == index[1]):
            behind_row = 1
            continue
        else:
            new_matrix[i - behind_row][replace] = 0.5 * (d_matrix[replace][i] + d_matrix[index[1]][i] - d_matrix[replace][index[1]])
            new_matrix[replace][i - behind_row] = new_matrix[i - behind_row][replace]
        
    return new_matrix


matrix_test = [[0, 5, 9, 9, 8], [5, 0, 10, 10, 9], [9, 10, 0, 8, 7], [9, 10, 8, 0, 3], [8, 9, 7, 3, 0]]

# gets the edges for the tree with neighbor joining algorithm
def get_tree(sequenceL_ist):
    distance_matrix = get_distances(sequenceL_ist)
    # first internal node made is number of tips + 1
    distance_id = len(distance_matrix) + 1

    # node ids list
    node_ids = []
    for i in range(1, len(distance_matrix) + 1):
        node_ids.append(i)

    # number of internal nodes is 
    # (number of iterations * 2) + 1
    edges = [[0 for _ in range(3)] for _ in range((len(distance_matrix)*2) - 3)]
    edge_index = 0
    final_iter = distance_id - 4

    # the iterations/steps of the algorithm are number of tips - 2 
    # that brings the matrix down to a 2x2 matrix from which the final joining can be done
    for i in range(0, distance_id - 3):
        # get the q matrix 
        
        x = calcQMatrix(distance_matrix)
        

        # get the index of first instance of minimum from q matrix
        # in format [x, y]

        z = findMin(x)

        # add joined nodes to edges matrix
        # first column has the node parent, second column has the node, third has distance between them

        b = distanceToNode(z, distance_matrix)

        edges[edge_index][0] = distance_id
        edges[edge_index][1] = node_ids[z[0]]
        edges[edge_index][2] = b[0]

        edge_index += 1

        edges[edge_index][0] = distance_id
        edges[edge_index][1] = node_ids[z[1]]
        edges[edge_index][2] = b[1]

        edge_index += 1

        # in node id list, replace smaller index with distance id to 
        # represent the new row/column
        # delete the larger index from node ids to represent row/column was deleted
        node_ids[z[0]] = distance_id
        node_ids.remove(node_ids[z[1]])

        if(i != final_iter):
            distance_id += 1
        # get new distance matrix with updated row and column/deleted row and column
        if(z[0] > z[1]):
          
            temp = z[0]
            z[0] = z[1]
            z[1] = temp
            
        # get new distance matrix with updated row and column/deleted row and column
        distance_matrix = updateDistanceMatrix(z, distance_matrix)

    final_dist = distance_matrix[0][1]

    # makes sure final root node is in correct item
    # since root node is the larger number(it's the last node to be made)
    if(node_ids[0] < node_ids[1]):
        temp = node_ids[1]
        node_ids[1] = node_ids[0]
        node_ids[0] = temp

    edges[edge_index][0] = node_ids[0]
    edges[edge_index][1] = node_ids[1]
    edges[edge_index][2] = final_dist

    return edges




# get the edges for the input fna file
edges = get_tree(sequence_list)


f.close()



f = open("edges.txt", "w")

# write the edges to edges.txt file
for x in edges:
    to_write = str(x[0]) + '\t' + str(x[1]) + '\t' + str(x[2])
    f.write(to_write)
    f.write('\n')

f.close()
# bootstrap support 

# this function takes in an edge matrix and returns the partitioning for
# each internal node
# format of nested dictionary, each node has it's own dictionary with left and right elem
# the left/right partitioning is represented as the addition of all ids for elements on left/right
# first encountered edge for node is left, second is right
# for example, take a node with edges [8, 4, 0.1] and [8, 5, 0.1]
# and a node with edges [9, 8, 0.3] and [9, 1, 0.2]
# partition would be {8: {left: 4, right: 5}, 9: {left: 9, right: 1}}


# original bootstrap implementation
# not used anymore but kept for reference


# def get_edges(edge_matrix):
#     partition_dict = {}
#     smallest = edge_matrix[0][0]
#     # left are even indices in edge list
#     # right are odd indices in edge list

#     # initialize the partition dictionary
#     for i in range(0, len(edge_matrix)):
#         if(edge_matrix[i][0] not in partition_dict):
#             partition_dict[edge_matrix[i][0]] = {"left": 0, "right": 0}
        
#         # special for root which has 3 possible children
#         # will not go out of bounds in check 
#         elif(edge_matrix[i][0] in partition_dict and i == len(edge_matrix) - 2):
#             partition_dict[edge_matrix[i][0]] = {"left": 0, "middle": 0, "right": 0}

#     for i in range(0, len(edge_matrix)):
#         # final node
#         if(i >= len(edge_matrix) - 3):
#             if(i == len(edge_matrix) - 3):
#                 partition_dict[edge_matrix[i][0]]["left"] = edge_matrix[i][1]
#                 if(edge_matrix[i][1] >= smallest):
#                     partition_dict[edge_matrix[i][0]]["left"] = partition_dict[edge_matrix[i][1]]["right"] + partition_dict[edge_matrix[i][1]]["left"]
#             elif(i == len(edge_matrix) - 2):
#                 partition_dict[edge_matrix[i][0]]["middle"] = edge_matrix[i][1]
#                 if(edge_matrix[i][1] >= smallest):
#                     partition_dict[edge_matrix[i][0]]["middle"] = partition_dict[edge_matrix[i][1]]["right"] + partition_dict[edge_matrix[i][1]]["left"]
#             else:
#                 partition_dict[edge_matrix[i][0]]["right"] = edge_matrix[i][1]
#                 if(edge_matrix[i][1] >= smallest):
#                     partition_dict[edge_matrix[i][0]]["right"] = partition_dict[edge_matrix[i][1]]["right"] + partition_dict[edge_matrix[i][1]]["left"]
#         elif(i % 2 == 0):
#             partition_dict[edge_matrix[i][0]]["left"] = edge_matrix[i][1]
#             if(edge_matrix[i][1] >= smallest):
#                     partition_dict[edge_matrix[i][0]]["left"] = partition_dict[edge_matrix[i][1]]["right"] + partition_dict[edge_matrix[i][1]]["left"]
        
#         else:
#             partition_dict[edge_matrix[i][0]]["right"] = edge_matrix[i][1]
#             if(edge_matrix[i][1] >= smallest):
#                     partition_dict[edge_matrix[i][0]]["right"] = partition_dict[edge_matrix[i][1]]["right"] + partition_dict[edge_matrix[i][1]]["left"]

    
#     return partition_dict




# this takes in an edge matrix and returns the same matrix
# but without the distances
# so [8, 5, 0.002] would become [8, 5] in this
def no_distances(matr):
    new_matr = []
    for x in matr:
        new_matr.append([x[0], x[1]])
    
    return new_matr


# this will get the tips associated with a node
# and put them in a list
# assumes that the given node is in the edge matrix
# also takes in the largest tip value to make sure it only returns
# tips and not nodes
# also takes in an empty list that will be added to
def get_descendant_tips(edge_matr, node, largest_tip, associated_tips):
    only_nodes = no_distances(edge_matr)
    for x in only_nodes:
        if(x[0] == node):
            if(x[1]<= largest_tip):
                associated_tips.append(x[1])
            # if it is an internal node, recursively get tips
            elif(x[1] > largest_tip):
                get_descendant_tips(edge_matr, x[1], largest_tip, associated_tips)
    
    return associated_tips



# this will get the partitioning of tips from a node and put them in a list
# for example if there is are edges like [[8, 5], [8, 4], [9, 1], [9, 8]]
# and partioning for node 8 is asked for
# it would return [[5], [4]]
# if partioning for node 9 is asked for
# it would return [[1], [5, 4]]
def get_partitioning(edge_matr, node, largest_tip, associated_partition):
    partition_set = set()
    only_nodes = no_distances(edge_matr)
    for x in only_nodes:
        associated_list = []
        if(x[0] == node):
            if(x[1] <= largest_tip):
                associated_partition.append([x[1]])
            else:
                # if internal node, get descendant tips for that node
                associated_partition.append(sorted(get_descendant_tips(edge_matr, x[1], largest_tip, associated_list)))

    # turn into sets because order for these doesn't matter
    # should make checking easier
    for x in associated_partition:
        set1 = set()
        for elem in x:
            set1.add(elem)
        
        partition_set.add(frozenset(set1))

    return partition_set



test_edges = [[22, 1], [22, 2], [23, 3], [23, 24], [24, 4], [24, 5], [21, 22], [21, 23]]
test_edges2 = [[38, 3], [38, 5], [37, 38], [37, 4], [36, 1], [36, 2], [35, 37], [35, 36]]
test_edges3 = [[19, 20], [19, 21], [20, 1], [20, 22], [21, 3], [21, 23], [22, 2], [22, 9], [23, 4], [23, 5]]
associated = []
associated2 = []
associated3 = []
x = get_partitioning(test_edges, 21, 10, associated)
y = get_partitioning(test_edges2, 35, 10, associated2)
z = get_partitioning(test_edges3, 19, 10, associated3)


def get_bootstraps(the_sequence_list):
    
    # this will be used to keep track of the largest tip
    largest = len(true_distance_matrix)


    # will use a dictionary of sets
    # basically each internal node in the original 
    # will have a dictionary entry with key being node id and value being the set of tips
    # then each bootstrap set of tips will be checked to see if there is a matching set in the dictionary
    # if so, that node will have its bootstrap confidence added to
    # also initializes the bootstrap confidence dictionary that will be added to

    found_nodes = []
    orig_dict = {}
    bootstrap_confidence = {}
    for x in edges:
        associated_list = []
        if(x[0] in found_nodes):
            continue
        else:
            orig_dict[x[0]] = (get_partitioning(edges, x[0], largest, associated_list))
            bootstrap_confidence[x[0]] = 0
            found_nodes.append(x[0])


    # will use to make bootstrap samples
    orig_sequence_list = copy.deepcopy(the_sequence_list)
    column_order = []
    
   
    # repeats for each bootstrap
    for i in range(0, int(bootstraps)):

        # replace each sequence column with random item from original 
        # chooses random number from (0, number of tips) so columns might repeat
        for l in range(0, len(the_sequence_list[0])):
            column_order.append(random.randint(0, len(the_sequence_list[0]) - 1))
        
        

        # for each string in sequence list
        # change to new column order

        # so for example if 3 is randomly chosen on the first column
        # it will replace all items from first column with all items from 3rd column
        # in original sequence list
        
        for j in range(0, len(the_sequence_list)):
            
            the_sequence_list[j] = list(the_sequence_list[j])
            
            for k in range(0, len(the_sequence_list[0])):
                the_sequence_list[j][k] = orig_sequence_list[j][column_order[k]]

            the_sequence_list[j] = ''.join(the_sequence_list[j])
        
        
        
        # columns have been arranged to make bootstrap samples
        # now make the tree with new sequence list
        # write edges file for each in a folder for making tree
        

        # file_name = 'bootstraps/bootstrap' + str(i) + ".txt"
        # f = open(file_name, "w")

        bootstrap_edges = get_tree(the_sequence_list)
        
       
        # for x in bootstrap_edges:
        #     to_write = str(x[0]) + '\t' + str(x[1]) + '\t' + str(x[2])
        #     f.write(to_write)
        #     f.write('\n')

        # f.close()

        # reset column order array 
        column_order = []
       

        
        # get partitioning for a certain node
        # check if that partitioning is present in the original tree
        # if so, increases confidence for the node that matched it
        found_nodes = []
        
        for x in bootstrap_edges:
            associated_list = []
            if x[0] in found_nodes:
                continue
            else:
                node_partition = get_partitioning(bootstrap_edges, x[0], largest, associated_list)
                for key in orig_dict:
                    if(node_partition == orig_dict[key]):
                        bootstrap_confidence[key] += 1
                        break
                found_nodes.append(x[0])
        
        # reset the sequence list
        the_sequence_list = copy.deepcopy(orig_sequence_list)
        
        
    return bootstrap_confidence
    

confidence = get_bootstraps(sequence_list)

f = open("bootstrap.txt", "w")

for key in confidence:
    to_write = str(key) + '\t' + str(confidence[key]/100)
    f.write(to_write)
    f.write('\n')
f.close()