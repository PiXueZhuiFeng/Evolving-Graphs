include("graph.jl")
using LinearAlgebra
using Random
using Statistics
using Laplacians
using SparseArrays
# using Threads
using DataFrames


mutable struct TreeNode
    parent :: Int                # 父节点,从那棵树更新过来的
    children::Vector{Int}        # 子节点的index
    node::Int                         # 节点，比上一个多加了一条边的边起点
    new_node::Int                     # 节点，比上一个多加了一条边的边终点
    weight :: Int                   
    root_rf :: Int              #rf的id        
end

mutable struct Rf
    next :: Vector{Int} # 每个点的父亲  长度n
    before :: Vector{Vector{Int}} ### 每个点的儿子们 n个list
    rootindex :: Vector{Int} #每个点的根节点 长度n
    children :: Vector{Int} ###更新之后可能的儿子们的index
    weight :: Int 
    id_thisforest :: Int ###这个森林的id
end





function add_edge_forest( i::Int, j::Int, allrf)
    ### useful_forest 存的是index
    ### newforest 存的是treenode的数组， 前面tau 是空的
    ### allrf 存的是所有的初始的森林
    for rf in allrf
        if useful_rf[rf.id_thisforest] == 1
            if  (rf.next[i] == -1 ) & (rf.rootindex[j]!=i)
                w = rf.weight
                new_child = TreeNode(rf.id_thisforest,Int[], i, j, w,rf.id_thisforest);   
                global id_newforest+=1; 
                push!(newforest,new_child)
                push!(rf.children,id_newforest);
                push!(useful_forest,id_newforest)
                ###原来的森林上加一个新的——end
            end
        end
    end
    ###原来的森林的儿子们上面加一个新的
    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        edge_list = Dict{Int, Int}()   
        while id_useful>tau
            # push!(edge_list,(current_forest.node,current_forest.new_node))
            edge_list[current_forest.node] = current_forest.new_node;
            id_useful = current_forest.parent;
            current_forest = newforest[id_useful];
        end
        ###结束之后id——useful就是初始的l个森林的id
        rf = allrf[id_useful]
        if (rf.next[i]==-1) & (get(edge_list,i,0)==0)   ###i符合要求，不是根
            root_j = rf.rootindex[j];
            while true
                newnode = get(edge_list,root_j,0);
                if newnode == 0
                    break
                else
                    root_j = rf.rootindex[newnode];
                end
            end
            if root_j !=i
                w = newforest[id_useful_i].weight
                new_child = TreeNode(id_useful_i,Int[], i, j, w,rf.id_thisforest);   
                global id_newforest+=1; 
                push!(newforest,new_child)
                push!(newforest[id_useful_i].children,id_newforest);
                push!(useful_forest,id_newforest)
            end
        end
    end
end



function del_edge_forest( i::Int, j::Int, allrf)
    ### useful_forest 存的是index
    ### newforest 存的是treenode的数组， 前面tau 是空的
    ### allrf 存的是所有的初始的森林
    ### useful_rf 是一个0,1数组，1表示这个森林还有用
    ### stat 是一个列表，他如果等于2表示这个森林可以被加上ij这条边，如果等于1表示他有ij这条边，如果等于3表示他不属于这两种情况，如果等于0表示他没有被查询到
    stat = zeros(Int,tau)
    for rf in allrf
        if useful_rf[rf.id_thisforest] == 1
            if  (rf.next[i] == -1 ) & (rf.rootindex[j]!=i)
                stat[rf.id_thisforest] =2;
                # rf.weight *= 1/2; 
                ###如果该森林可以被加上ij这条边，那么就给他的权重不变       
            elseif rf.next[i] == j 
                ###如果该森林有ij这条边，那么就给他的权重不变. 删掉ij这条边，并更新他的before
                stat[rf.id_thisforest] =1;
                rf.next[i] = -1;
                idx = findfirst(==(i), rf.before[j])
                deleteat!(rf.before[j], idx)
                stack = [i]  # 初始化栈，起始包含根节点
                
                while !isempty(stack)
                    node = pop!(stack)  # 取出当前节点
                    rf.rootindex[node] = i; #改他的rootindex
                    # 将子节点加入栈中
                    for child in (rf.before[node])
                        push!(stack, child)
                    end
                end
                # rf.weight *= 1/2; 
            else
                ### 剩下的情况，权重乘2
                 rf.weight *= 2; 
                stat[rf.id_thisforest] =3;
            end
        end
    end
    a1 = 0
    a2 = 0
    a3 = 0
    for i in stat
        if i == 1
            a1+=1;
        elseif i==2
            a2+=1;
        elseif i==3
            a3+=1
        end
    end

    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        rf_index = current_forest.root_rf;
        if stat[rf_index] == 0
            rf = allrf[rf_index]
            if  (rf.next[i] == -1 ) & (rf.rootindex[j]!=i)
                stat[rf.id_thisforest] =2;
                ###如果该森林可以被加上ij这条边，那么就给他的权重不变     
            elseif rf.next[i] == j 
                ###如果该森林有ij这条边，那么就给他的权重不变. 删掉ij这条边，并更新他的before
                stat[rf.id_thisforest] =1;
                rf.next[i] = -1;
                idx = findfirst(==(i), rf.before[j])
                deleteat!(rf.before[j], idx)
                stack = [i]  # 初始化栈，起始包含根节点
                
                while !isempty(stack)
                    node = pop!(stack)  # 取出当前节点
                    rf.rootindex[node] = i; #改他的rootindex
                    # 将子节点加入栈中
                    for child in (rf.before[node])
                        push!(stack, child)
                    end
                end
            else
                ### 剩下的情况，权重乘2
                rf.weight *= 2;
                stat[rf.id_thisforest] =3;
            end
        end


        if stat[rf_index] == 1  #原来的森林有ij这条边并且已经删掉了，这时候不需要操作
        elseif stat[rf_index] == 2 #原来的森林可以加ij这条边，这时候需要再check一下
            id_useful = id_useful_i
            edge_list = Dict{Int, Int}()   
            while id_useful>tau
                edge_list[current_forest.node] = current_forest.new_node;
                id_useful = current_forest.parent;
                current_forest = newforest[id_useful];
            end
            ###下面找ij的根
            rf = allrf[id_useful]
            root_j = rf.rootindex[j];
            while true
                newnode = get(edge_list,root_j,0);
                if newnode == 0
                    break
                else
                    root_j = rf.rootindex[newnode];
                end
            end

            root_i = rf.rootindex[i];
            while true
                newnode = get(edge_list,root_i,0);
                if newnode == 0
                    break
                else
                    root_i = rf.rootindex[newnode];
                end
            end
            ####寻找结束
            id_useful = id_useful_i
            current_forest = newforest[id_useful];
            if (i==root_i) & (root_j!=i)
            elseif (get(edge_list,i,0)==j)
                id_useful = id_useful_i
                current_forest = newforest[id_useful];
                while !((current_forest.node == i) & (current_forest.new_node == j))
                    id_useful = current_forest.parent;
                    current_forest = newforest[id_useful];
                end
                current_forest.node = G.n+1;
                current_forest.new_node = G.n+1;
            else
                current_forest.weight  *=2;
            end
        else
            id_useful = id_useful_i
            current_forest = newforest[id_useful];
            current_forest.weight  *=2;
        end
    end
end


function searchij(i,j,allrf, newforest, useful_forest)
####查i的根是j的概率
    ansij = 0; 
    sumw = 0;
    for rf in allrf 
        if useful_rf[rf.id_thisforest] == 1
            wrf =  rf.weight
            sumw += wrf ;
            if rf.rootindex[i] == j 
                ansij += wrf;
            end
        end
    end
    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        edge_list = Dict{Int, Int}()   
        while id_useful>tau
            # push!(edge_list,(current_forest.node,current_forest.new_node))
            edge_list[current_forest.node] = current_forest.new_node;
            id_useful = current_forest.parent;
            current_forest = newforest[id_useful];
        end
        ###结束之后id——useful就是初始的l个森林的id
        rf = allrf[id_useful]
        # if (rf.next[i]==-1) & (get(edge_list,i,0)==0)   ###i符合要求，不是根
            root_i = rf.rootindex[i];
            while true
                newnode = get(edge_list,root_i,0);
                if newnode == 0
                    break
                else
                    root_i = rf.rootindex[newnode];
                end
            end
            wrf = newforest[id_useful_i].weight
            sumw += wrf ;
            if root_i ==j
                ansij += wrf;
            end
        # end
    end
    return ansij/sumw
end


function searchijplus(i,j,allrf, newforest, useful_forest,d,L)
    ####查i的根是j的概率
    sumw = 0;
    ansij = 0;
    for rf in allrf 
        if useful_rf[rf.id_thisforest] == 1
            wrf =  rf.weight
            sumw += wrf ;
            k  = rf.rootindex[i];
            if (k == j) & (i != j )
                ansij += wrf/(2+d[j])
            elseif L[k,j] <0 
                if i== j
                    ansij += wrf/(1+d[j])
                else
                    ansij += wrf/(2+d[j])
                end
            end
        end
    end
    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        edge_list = Dict{Int, Int}()   
        while id_useful>tau
            edge_list[current_forest.node] = current_forest.new_node;
            id_useful = current_forest.parent;
            current_forest = newforest[id_useful];
        end
        ###结束之后id——useful就是初始的l个森林的id
        rf = allrf[id_useful]
        root_i = rf.rootindex[i];
        while true
            newnode = get(edge_list,root_i,0);
            if newnode == 0
                break
            else
                root_i = rf.rootindex[newnode];
            end
        end
        wrf = newforest[id_useful_i].weight
        k = root_i;
        sumw += wrf ;
        if (k == j) & (i != j )
            ansij += wrf/(2+d[j])
        elseif L[k,j] <0 
            if i== j
                ansij += wrf/(1+d[j])
            else
                ansij += wrf/(2+d[j])
            end
        end
    end
    ansij = ansij/sumw 
    if i == j 
       ansij += 1/(1+d[i]);
    end
    return ansij
end    







function searchirow(i,allrf, newforest, useful_forest)
    sumw = 0;
    id_j = Int[];
    val_j = Float32[];
    for rf in allrf 
        if useful_rf[rf.id_thisforest] == 1
        wrf = rf.weight
        sumw += wrf ;
        push!(id_j,rf.rootindex[i])
        push!(val_j,wrf)
        end
    end
    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        edge_list = Dict{Int, Int}()   
        while id_useful>tau
            edge_list[current_forest.node] = current_forest.new_node;
            id_useful = current_forest.parent;
            current_forest = newforest[id_useful];
        end
        ###结束之后id——useful就是初始的l个森林的id
        rf = allrf[id_useful]
        root_i = rf.rootindex[i];
        while true
            newnode = get(edge_list,root_i,0);
            if newnode == 0
                break
            else
                root_i = rf.rootindex[newnode];
            end
        end
        w = newforest[id_useful_i].weight
        wrf = w
        sumw += wrf ;
        push!(id_j,root_i)
        push!(val_j,wrf)
    end
    val_j = val_j./sumw
    row_i = sparsevec(id_j,val_j,G.n)
    return row_i
end

function searchirowplus(i,allrf, newforest, useful_forest,d)
    ####查第i行
    sumw = 0;
    id_j = Int[];
    val_j = Float32[];
    for rf in allrf 
        if useful_rf[rf.id_thisforest] == 1
            wrf =  rf.weight 
            sumw += wrf ;
            k = rf.rootindex[i]
            nn = nbr[k];
            for s in nn 
                if i == s
                    push!(id_j,i)
                    push!(val_j,wrf/(1+d[i]))
                else
                    push!(id_j,s) 
                    push!(val_j,wrf/(2+d[s]))
                end
            end
            if i != k
                push!(id_j,k)
                push!(val_j,wrf/(2+d[k]))
            end
        end
    end
    for id_useful_i in useful_forest
        id_useful = id_useful_i
        current_forest = newforest[id_useful];
        edge_list = Dict{Int, Int}()   
        while id_useful>tau
            edge_list[current_forest.node] = current_forest.new_node;
            id_useful = current_forest.parent;
            current_forest = newforest[id_useful];
        end
        ###结束之后id——useful就是初始的l个森林的id
        rf = allrf[id_useful]
        root_i = rf.rootindex[i];
        while true
            newnode = get(edge_list,root_i,0);
            if newnode == 0
                break
            else
                root_i = rf.rootindex[newnode];
            end
        end
        w = newforest[id_useful_i].weight
        wrf =  w
        sumw += wrf ;

        k = root_i
        nn = nbr[k];
        for s in nn 
            if i == s
                push!(id_j,i)
                push!(val_j,wrf/(1+d[i]))
            else
                push!(id_j,s) 
                push!(val_j,wrf/(2+d[s]))
            end
        end
        if i != k
            push!(id_j,k)
            push!(val_j,wrf/(2+d[k]))
        end
    end
    val_j = val_j./sumw;
    push!(id_j,i)
    push!(val_j,1/(1+d[i]))
    row_i = sparsevec(id_j,val_j,G.n)
    return row_i
end




 
 

 



function rfsave(G,tau)
    t1 = time()
    #### 初始化
    d = [length(G.nbr[j]) for j = 1:G.n];
    InForest = [false for j = 1:G.n];
    # allrf = [Rf(Int64[], [Int64[]], Int64[], Int64[], 0, 0) for i = 1:tau];
    # allrf = Rf[];
    GC.enable(false)
    allrf = [Rf(-ones(UInt32,G.n),   [UInt32[] for j=1:G.n+1], zeros(UInt32, G.n), TreeNode[], 1, 0) for i = 1:tau];
    GC.enable(true)
    GC.gc()
    # Next = -ones(Int,G.n);
    # rootindex = zeros(Int, G.n);
    GC.enable(false)
    for p = 1:tau
        InForest .= false;
        Next = allrf[p].next;
        rootindex = allrf[p].rootindex;
        # Next .=-1
        # rootindex.=0
        for i = 1:G.n
            u = i;
            while !InForest[u]
                seed = rand();
                if seed <= 1 / (1 + d[u])
                    InForest[u] = true;
                    Next[u] = -1;
                    rootindex[u] = u;
                else
                    k = floor(Int, seed * (1 + d[u]));
                    Next[u] = G.nbr[u][k];
                    u = Next[u];   
                end
            end
            rootnow = rootindex[u];
            u = i;
            while !InForest[u]
                InForest[u] = true;
                rootindex[u] = rootnow;
                u = Next[u];
            end
        end
        bef = allrf[p].before;
        for i=1:G.n
            j = Next[i]
            if j== -1
                push!(bef[G.n+1],i)
            else
                push!(bef[j],i)
            end
        end
        allrf[p].id_thisforest = p
        # rf = Rf(Next,bef,rootindex,TreeNode[],1,p);
        # push!(allrf,rf)
        # allrf[p] = rf
    end
    GC.enable(true)
    GC.gc()
    t2 = time()
    println(t2-t1)
    return allrf
end

function rfsavemul(G,tau)
    t1 = time()
    #### 初始化
    d = [length(G.nbr[j]) for j = 1:G.n];
    # InForest = [false for j = 1:G.n];
    # allrf = [Rf(Int64[], [Int64[]], Int64[], Int64[], 0, 0) for i = 1:tau];
    # allrf = Rf[];
    allrf = Vector{Rf}(undef, tau)
    GC.enable(false)
    Threads.@threads for i = 1:tau
        allrf[i] = Rf(-ones(UInt32, G.n), [UInt32[] for j=1:G.n+1], zeros(UInt32, G.n), TreeNode[], 1, 0)
    end
    GC.enable(true)
    GC.gc()
    # Next = -ones(Int,G.n);
    # rootindex = zeros(Int, G.n);
    GC.enable(false)
    Threads.@threads for p = 1:tau
        InForest = [false for j = 1:G.n];
        Next = allrf[p].next;
        rootindex = allrf[p].rootindex;
        # Next .=-1
        # rootindex.=0
        for i = 1:G.n
            u = i;
            while !InForest[u]
                seed = rand();
                if seed <= 1 / (1 + d[u])
                    InForest[u] = true;
                    Next[u] = -1;
                    rootindex[u] = u;
                else
                    k = floor(Int, seed * (1 + d[u]));
                    Next[u] = G.nbr[u][k];
                    u = Next[u];   
                end
            end
            rootnow = rootindex[u];
            u = i;
            while !InForest[u]
                InForest[u] = true;
                rootindex[u] = rootnow;
                u = Next[u];
            end
        end
        bef = allrf[p].before;
        for i=1:G.n
            j = Next[i]
            if j== -1
                push!(bef[G.n+1],i)
            else
                push!(bef[j],i)
            end
        end
        allrf[p].id_thisforest = p
        # rf = Rf(Next,bef,rootindex,TreeNode[],1,p);
        # push!(allrf,rf)
        # allrf[p] = rf
    end
    GC.enable(true)
    GC.gc()
    t2 = time()
    println(t2-t1)
    return allrf
end



 
function cutbranch(tt,useful_forest,useful_rf)
    tot = 0;
    for rf in allrf
        if useful_rf[rf.id_thisforest] ==1
            tot+= rf.weight
        end
    end
    for i in useful_forest
        tot+=   newforest[i].weight
    end
    tot /=tt;
    for rf in allrf
        if useful_rf[rf.id_thisforest] ==1
            w = rf.weight
            neww = 0;
            for j = 1:w 
                if rand()<1 /tot
                    neww += 1;
               end
            end
            rf.weight = neww;
            if neww == 0
                useful_rf[rf.id_thisforest] =0
            end
        end
    end
    new_useful_forest = Int[];
    for i in useful_forest
        w = newforest[i].weight
        neww = 0;
        for j = 1:w 
            if rand()<1 /tot
                neww += 1;
            end
        end
        newforest[i].weight = neww;
        if neww != 0
            push!(new_useful_forest,i)
        end
    end
    global useful_forest = new_useful_forest;
end



 

 