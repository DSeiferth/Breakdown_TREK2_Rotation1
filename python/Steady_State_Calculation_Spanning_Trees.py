import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
#import Bib3_version4b as version4b
import time



def _expand(G, explored_nodes, explored_edges):
    """
    Expand existing solution by a process akin to BFS.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    explored_nodes: set of ints
        nodes visited

    explored_edges: set of 2-tuples
        edges visited

    Returns:
    --------
    solutions: list, where each entry in turns contains two sets corresponding to explored_nodes and explored_edges
        all possible expansions of explored_nodes and explored_edges

    """
    frontier_nodes = list()
    frontier_edges = list()
    for v in explored_nodes:
        for u in nx.neighbors(G,v):
            if not (u in explored_nodes):
                frontier_nodes.append(u)
                frontier_edges.append([(u,v), (v,u)])

    return zip([explored_nodes | frozenset([v]) for v in frontier_nodes], [explored_edges | frozenset(e) for e in frontier_edges])

def find_all_spanning_trees(G, root=0):
    """
    Find all spanning trees of a Graph.

    Arguments:
    ----------
    G: networkx.Graph() instance
        full graph

    Returns:
    ST: list of networkx.Graph() instances
        list of all spanning trees

    """

    # initialise solution
    explored_nodes = frozenset([root])
    explored_edges = frozenset([])
    solutions = [(explored_nodes, explored_edges)]
    # we need to expand solutions number_of_nodes-1 times
    for ii in range(G.number_of_nodes()-1):
        # get all new solutions
        solutions = [_expand(G, nodes, edges) for (nodes, edges) in solutions]
        # flatten nested structure and get unique expansions
        solutions = set([item for sublist in solutions for item in sublist])

    return [nx.from_edgelist(edges) for (nodes, edges) in solutions]


def Create_Graph(N, seed, ausgabe):
    H=nx.grid_2d_graph(N, 2, periodic=False, create_using=nx.DiGraph())
    #print(list(H.nodes))

    labels = dict( ((i,j), 2*i+j) for i, j in H.nodes() )
    #print labels

    nx.relabel_nodes(H,labels,False)
    inds=labels.keys()
    vals=labels.values()
    inds=[(N-j-1,N-i-1) for i,j in inds]
    pos2=dict(zip(vals,inds))

    fig, ax = plt.subplots(1,1)
    nx.draw_networkx(H, pos=pos2, with_labels=True, node_size = 200, node_color='orange',font_size=10,ax=ax)
    plt.axis('off')
    plt.title('grid')
    plt.gca().invert_xaxis()
    if ausgabe==0:
        plt.show()
    
    random.seed(seed)
    for u,v,d in H.edges(data=True):
        H[u][v]['weight']=random.uniform(0.0, 4.0)
        #d['weight']+=7
    if ausgabe==0:
        print(H.edges.data('weight'))
    return [H, pos2]

def weights(G,s,t):
    N=G.number_of_nodes()
    #for nodes in G:
    bla=[] 
    for path in nx.all_simple_paths(G, source=s, target=t):
        bla.append(path)
    print(bla)
    weight=1.0
    pfad=bla[0]
    for i in range(len(pfad)-1):
        print(G[pfad[i]][pfad[i+1]]['weight'])
        weight=weight*G[pfad[i]][pfad[i+1]]['weight']
    return  weight
    
def calc_weights(G,ZIEL):
    ausgabe=1
    N=G.number_of_nodes()
    # Erstelle Liste mit Pfaden, nehme den laengsten
    Pfadliste=[]
    for node in G:
        if node!=ZIEL:
            bla=[]
            for path in nx.all_simple_paths(G, source=node, target=ZIEL):
                bla.append(path)
            if len(bla)>1:
                print('mehr als einen path in calc_weights')
                print(bla)
            Pfadliste.append(bla[0])
    #print Pfadliste
    
    #Finde die laengeste Liste
    length=np.array([])
    for i in Pfadliste:
        length=np.append(length, len(i))
    pos = np.where(length == np.amax(length))[0][0]
    Pfad = Pfadliste[pos]
    
    #Berechne Gewicht des Pfads
    weight=1.0
    for i in range(len(Pfad)-1):
        #print G[Pfad[i]][Pfad[i+1]]['weight']
        weight=weight*G[Pfad[i]][Pfad[i+1]]['weight']
    
    # Finde nodes, die nicht auf Pfad liegen
    missing=[]
    for node in G:
        if node not in Pfad:
            #print node
            missing.append(node)
  
    # fehlende Pfade in einer Liste
    missingpaths=[]
    for m in missing:
        bla=[]
        for path in nx.all_simple_paths(G, source=m, target=ZIEL):
            bla.append(path)
        missingpaths.append( bla[0])
    
        
    calculatedPaths=[Pfad]
    
    # Gewichte der fehlende Pfade berechnen
    # dabei aber keine edge doppelt beruecksichtigen
    faktor=1.0
    while(len(missingpaths)>0):
        #finde laengstes Element
        length=np.array([])
        for i in missingpaths:
            length=np.append(length, len(i))
        pos = np.where(length == np.amax(length))[0][0]
        pfad2 = missingpaths[pos]   # pfad2 ist der laengste pfad in misingpaths
        
        i=0
        for m in pfad2:
            use=0
            for schonDa in calculatedPaths:
                if m in schonDa:
                    use=1
            if use==0:
                faktor=faktor*G[pfad2[i]][pfad2[i+1]]['weight']
            i=i+1
        
        missingpaths.remove(pfad2)  # laengster pfad wird entfernt
        calculatedPaths.append(pfad2)
        # sind andere pfade bereits enthalten
        for bla in missingpaths:
            if bla[0] in pfad2:
                missingpaths.remove(bla)
                
    return faktor*weight
        
    

def weight_of_tree(G):
    #gibt Gewicht fuer einen Baum fuer jeden zustand
    weight=[]
    #for node in G:
    for node in range(len(G)):
        #print('weight of tree: directed towards node '+str(node))
        #print(node)
        fak=calc_weights(G,node)
        weight.append(fak)
        #print 
    return weight
        
    
def steady_state(G):
    Gewichte=[] 
    ST = find_all_spanning_trees(G)
    #print('Maximal trees')
    for g in ST:
        #print (g.edges.data())
        g=g.to_directed()
        for node1 in g:
            for node2 in g:
                if g.has_edge(node1, node2):
                    g[node1][node2]['weight']=G[node1][node2]['weight']
        #print (g.edges.data())
        Gewichte.append(weight_of_tree(g) )
    #print('Gewichtsvektor für maximales trees')
    #print(Gewichte)
    #Wahrscheinlichkeiten ausrechnen
    NoRows=len(Gewichte)
    NoColums=len(Gewichte[0])
    p=[]
    for i in range(NoColums):
        p_i=0
        for tree in range(NoRows):
            p_i=p_i+Gewichte[tree][i]
        p.append(p_i)
    #print p
    Normfactor=0.0
    for q in p:
        Normfactor=Normfactor+q
    p=np.array(p)*1.0/Normfactor
    #print p 
    return p

def Normfactor(G):
    Gewichte=[]
    ST = find_all_spanning_trees(G)
    for g in ST:
        g=g.to_directed()
        for node1 in g:
            for node2 in g:
                if g.has_edge(node1, node2):
                    g[node1][node2]['weight']=G[node1][node2]['weight']
        #print g.edges.data()
        Gewichte.append(weight_of_tree(g) )
    #print Gewichte
    #Wahrscheinlichkeiten ausrechnen
    NoRows=len(Gewichte)
    NoColums=len(Gewichte[0])
    p=[]
    for i in range(NoColums):
        p_i=0
        for tree in range(NoRows):
            p_i=p_i+Gewichte[tree][i]
        p.append(p_i)
    #print p
    Normfactor=0.0
    for q in p:
        Normfactor=Normfactor+q
    p=np.array(p)*1.0/Normfactor
    return Normfactor


def Graph2Matrix(G):
    N=nx.number_of_nodes(G)
    A=np.zeros((N, N))
    for node1 in G:
        for node2 in G:
            if G.has_edge(node1, node2):
                A[node1][node2]=G[node1][node2]['weight']
    #print A
    return A

def Matrix2Graph(a):
    N=len(a)
    G = nx.DiGraph()
    for row in range(N):
        for col in range(N):
            if a[row][col]!=0:
                G.add_weighted_edges_from([(row, col, a[row][col] )])
    #print G.edges(data=True)
    return G

def simulation_prob(matrix, startState, acc):
    start_time = time.time()
    [sumUp, rateInterval]=version4b.init(matrix)
    t=0.0
    state=startState
    #memory for trajectory
    TIME=[0]
    states=[startState]
    state_time=np.zeros(len(matrix))
    prob=np.zeros(len(matrix))
    prob_prev=np.zeros(len(matrix))
    i=0
    diff=1.0
    while diff>acc:
        [state, t]=version4b.timestep(state, t, sumUp, rateInterval)
        #print state
        #print t
        TIME.append(t)
        states.append(state)
        dt=TIME[i+1]-TIME[i]
        state_time[ states[i] ]=state_time[ states[i] ]+dt
        
        prob_prev=prob
        prob=1.0/t*state_time
        #print prob_prev
        #print prob
        i=i+1
        diff=max(abs( np.array(prob) - np.array(prob_prev) ))
    print("--- %s seconds ---for simulation prob" % (time.time() - start_time))
    return prob
    
def control(G, acc):
    print('Analytical steady state')
    print(steady_state(G))
    print
    matrix=Graph2Matrix(G)
    print(matrix)
    print('Simulated steady state')
    #print version4b.simulation_prob(matrix, 0, T)
    print(simulation_prob(matrix, 0, acc))
