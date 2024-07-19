#!/usr/bin/env python3

from my_imports import *

def table(k):
    mytable=[[0],[1]]
    if k>1:
        for i in range(k-1):
            # double_table=[]
            add_ref=[]
            for list in mytable:
                newlist=list[:]
                newlist.append(0)
                add_ref.append(newlist)
            add_alt=[]
            for list in mytable:
                newlist=list[:]
                newlist.append(1)
                add_alt.append(newlist)  
            double_table=add_ref+add_alt
            mytable=double_table
        mytable.remove([0]*k)
        mytable.remove([1]*k)
    return mytable

def init_delta(points_num):
    delta_set=[]
    for i in range(points_num-1):
        delta_set.append([[0,0],[0,0]])
    return delta_set

def fixed(k):
    num_list=[]
    for i in range(k):
        num_list.append(i+1)
    list_sum=sum(num_list)
    fixed_alpha=[]
    for num in num_list:
        fixed_alpha.append(num/list_sum)
    return fixed_alpha
    # alpha = np.array([np.random.randint(1, 10) for i in range(k)])
    # alpha = np.sort(alpha)
    # fixed_alpha= (alpha/sum(alpha)).tolist()
    # return fixed_alpha

def alpha_step(delta_set,geno_set,k,answer_index,beta_set,share_set,weight): 
    prob = LpProblem('myPro', LpMinimize)
    alpha=LpVariable.dicts("alpha",range(k),lowBound = 0,upBound = 1)
    cost=LpVariable.dicts("cost",range(len(geno_set)*4),lowBound = 0)
    beta_cost=LpVariable.dicts("beta_cost",range(len(geno_set)),lowBound = 0)
    alpha_sum=0
    for i in range(k):
        alpha_sum+=alpha[i]
    prob+=alpha_sum==1
    #cost for delta
    delta_lost=0
    beta_lost=0
    for i in range(len(delta_set)):

        supposed_delta=[[0,0],[0,0]]
        for l in range(k):
            supposed_delta[geno_set[i][l]][geno_set[i+1][l]] += alpha[l]
        prob+=delta_set[i][0][0] - supposed_delta[0][0] <= cost[4*i]
        prob+=delta_set[i][0][0] - supposed_delta[0][0] >=-cost[4*i]
        delta_lost+=cost[4*i]*(1-weight)
        prob+=delta_set[i][0][1] - supposed_delta[0][1] <= cost[4*i+1]
        prob+=delta_set[i][0][1] - supposed_delta[0][1] >=-cost[4*i+1]
        delta_lost+=cost[4*i+1]*(1-weight)
        prob+=delta_set[i][1][0] - supposed_delta[1][0] <= cost[4*i+2]
        prob+=delta_set[i][1][0] - supposed_delta[1][0] >=-cost[4*i+2]
        delta_lost+=cost[4*i+2]*(1-weight)
        prob+=delta_set[i][1][1] - supposed_delta[1][1] <= cost[4*i+3]
        prob+=delta_set[i][1][1] - supposed_delta[1][1] >=-cost[4*i+3]
        delta_lost+=cost[4*i+3]*(1-weight)
    #cost for beta
    for j in range(len(beta_set)):
        myratio=0
        for l in range(k):
            myratio+=geno_set[j][l]*alpha[l]
        # beta=delta_set[i][0][1] + delta_set[i][1][1]
        beta=beta_set[j]
        prob+=myratio-beta <=beta_cost[j]
        prob+=myratio-beta >=-beta_cost[j]
        beta_lost+=beta_cost[j]*weight

    lost=0
    lost+=delta_lost
    lost+=beta_lost

    prob+=lost,"total lost"
    prob.solve()

    # alpha = []
    # for i in prob.variables():
    #     if 'alpha' in i.name:
    #          alpha.append(i.varValue)

    alpha = np.empty(k) # avoid index ambiguity
    for i in prob.variables():
        if 'alpha' in i.name:
            index = int(i.name.split("_")[1])
            # print (i.name, i.varValue)
            alpha[index] = i.varValue
    # print (alpha)

    return alpha

def index2seq(strain_number,locus_index):
    allele_set=[]
    my_table=table(strain_number)
    # print (my_table,locus_index)
    for index in locus_index:
        allele_set.append(my_table[index])
    allele_set=np.array(allele_set)
    seq_list=np.transpose(allele_set)
    return seq_list

class Phase_step(): #too many same step in the iteration
    def __init__(self,given_alpha,delta_set,beta_set,share_set,weight):
        self.beta_set=beta_set
        self.k=len(given_alpha)
        self.given_alpha=given_alpha
        self.delta_set=delta_set
        self.points_num=len(delta_set)+1
        self.table=table(self.k)
        self.table_delta=self.delta_table()
        self.estimated_beta=self.estimate()
        # self.delta_exp=self.expected_delta()
        self.share_set=share_set
        self.w=weight
    def delta_diff(self,delta_a,delta_b):
        return abs(delta_a[0][0]-delta_b[0][0])+abs(delta_a[0][1]-delta_b[0][1])+\
            abs(delta_a[1][0]-delta_b[1][0])+abs(delta_a[1][1]-delta_b[1][1])
    def delta_table(self):
        num=len(self.table)
        table_delta=[]
        for i in range(num):
            middle_table=[]
            for j in range(num):
                ij_delta=[[0.0,0.0],[0.0,0.0]]
                for e in range(self.k):
                    ij_delta[self.table[i][e]][self.table[j][e]]+=self.given_alpha[e]
                middle_table.append(ij_delta)
            table_delta.append(middle_table)
        return table_delta
    def delta_phase(self):
        # self.delta_set[start:end+1],self.beta_set[start:end+2]=graph.lp_correct(self.delta_set[start:end+1],self.beta_set[start:end+2])
        save_table=[]
        geno_num=len(self.table)
        for r in range(len(self.delta_set)+1): 
            point_table=[]
            # print ('r',r,self.share_set[r-1],2**(int(-self.share_set[r-1])))
            if r == 0:
                for m in range(geno_num):
                    beta_loss=abs(self.estimated_beta[m]-self.beta_set[r])
                    weight_loss=beta_loss*(self.w)
                    point_table.append([weight_loss,0])
                save_table.append(point_table)
            else:
                for m in range(geno_num):
                    this_geno=[float('inf'),0]
                    for n in range(geno_num): 
                        # print ('sadfsfd', self.delta_set[r-1])
                        beta_loss=abs(self.estimated_beta[m]-self.beta_set[r])
                        # if sum(self.delta_set[r-1][0]) + sum(self.delta_set[r-1][1]) == 0:
                        #     delta_loss = beta_loss*(1 - self.w)
                        # else:
                        #     delta_loss=self.delta_diff(self.delta_set[r-1],self.table_delta[n][m])
                        delta_loss=self.delta_diff(self.delta_set[r-1],self.table_delta[n][m])
                        
                        weight_loss=beta_loss*(self.w)+delta_loss*(1-self.w)
                        update_loss=save_table[r-1][n][0] + weight_loss
                        add_loss=round(update_loss,6)
                        if add_loss<this_geno[0]:
                            this_geno=[add_loss,n]
                    point_table.append(this_geno)
                save_table.append(point_table)
        frag_index,part_loss=self.backtrack(save_table)
        # print ('save table', save_table, part_loss, frag_index, self.table)
        # delta_sum=self.loss(frag_index)
        # return frag_index,delta_sum
        return frag_index,part_loss       
    def backtrack(self,save_table):
        # geno_num=2**self.k-2
        geno_num=len(self.table)
        reverse_index=[]
        final_geno=[float('inf'),0]
        final_index=0
        for m in range(geno_num):
            this_geno=save_table[-1][m]
            if float(this_geno[0]) < float(final_geno[0]):
                final_geno=this_geno
                final_index=m
        # print ("delta loss",final_geno[0])
        part_loss=final_geno[0]
        reverse_index.append(final_index)
        for r in reversed(range(len(save_table)-1)):
            reverse_index.append(final_geno[1])
            final_geno=save_table[r][reverse_index[-1]]
        reverse_index.reverse()
        answer_index=reverse_index
        return answer_index,part_loss 
    def breaks_phase(self):
        answer_index,phase_loss=self.delta_phase()
        geno_set=self.genotype(answer_index)
        return answer_index,geno_set,phase_loss
    def genotype(self,answer_index):
        geno_set=[]
        
        for i in range(self.points_num):
            geno_set.append(self.table[int(answer_index[i])]) 
        return geno_set
    def estimate(self):
        estimated_beta=[]
        for colum in self.table:
            result=sum(np.multiply(np.array(colum),np.array(self.given_alpha)))
            estimated_beta.append(round(result,6))
        return estimated_beta
    def expected_delta(self):
        delta_exp=[]
        for i in range(len(self.estimated_beta)):
            middle_table=[]
            for j in range(len(self.estimated_beta)):
                fir=self.estimated_beta[i]
                sec=self.estimated_beta[j]
                exp=[[(1-fir)*(1-sec),(1-fir)*sec],[fir*(1-sec),fir*sec]]
                middle_table.append(exp)
            delta_exp.append(middle_table)
        return delta_exp
    def loss(self,answer_index):
        delta_sum=0
        for r in range(len(self.delta_set)): 
            delta_loss=self.delta_diff(self.delta_set[r],self.table_delta[answer_index[r]][answer_index[r+1]])
            # delta_sum+=delta_loss*(1-2**(int(-self.share_set[r])))
            if int(self.share_set[r])>4:
                delta_sum+=delta_loss
        return delta_sum

class Workflow():
    def __init__(self,beta_set,delta_set,share_set,weight,elbow):
        # print (delta_set)
        self.delta_set=delta_set
        self.beta_set=beta_set
        self.share_set=share_set
        self.w=weight
        self.elbow=elbow
    def choose_k(self):
        previous_loss=float('inf')
        T=1
        while True:
            geno_index, corr_loss, final_alpha=self.iteration(T)
            #print (T, corr_loss, corr_loss/len(self.beta_set) , float(previous_loss-corr_loss)/previous_loss)
            if corr_loss==0 or float(previous_loss-corr_loss)/previous_loss < float(self.elbow) or corr_loss/len(self.beta_set) < 0.02:
                if T == 1:
                    previous_alpha,previous_index=final_alpha,geno_index
                elif float(previous_loss-corr_loss)/previous_loss >= float(self.elbow):
                    previous_alpha,previous_index=final_alpha,geno_index
                seq_list=index2seq(len(previous_alpha),previous_index)
                return previous_alpha,seq_list,corr_loss
            previous_loss=corr_loss
            previous_alpha,previous_index=final_alpha,geno_index
            T+=1        
    def given_k(self,T):
        geno_index,corr_loss,final_alpha=self.iteration(T)
        seq_list=index2seq(len(final_alpha),geno_index)
        return final_alpha,seq_list,corr_loss
    def iteration(self,T):
        #T is supposed strain number
        times=0
        past_loss=float('inf')
        save_list=[]
        loss_list=[]
        current_alpha=fixed(T)    
        while True:  
            ph=Phase_step(current_alpha,self.delta_set,self.beta_set,self.share_set,self.w)
            answer_index,geno_set,phase_loss=ph.breaks_phase()
            save_list.append([answer_index,geno_set,phase_loss,current_alpha])
            loss_list.append(phase_loss)
            if phase_loss==0 or abs(past_loss-phase_loss) <0.000001 or times > 15:
                final_index,final_loss,final_alpha = answer_index,phase_loss,current_alpha
                break
            past_loss=phase_loss
            times+=1
            if T==1:
                current_alpha=[1.0]
            else:
                current_alpha=alpha_step(self.delta_set,geno_set,T,answer_index,self.beta_set,self.share_set,self.w)
            current_alpha=sorted(current_alpha)
        return final_index,final_loss,final_alpha
