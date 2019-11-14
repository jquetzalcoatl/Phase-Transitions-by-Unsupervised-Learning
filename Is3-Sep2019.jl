using LinearAlgebra, DelimitedFiles, Arpack, IterativeSolvers,PlotlyJS
#plotlyjs() Plots
mutable struct Site
    snum
    spin
    vec
    z
end
Ldim=40;
J=-1;
T=10;
snap=100;
sites=Array{Site}(undef,Ldim,Ldim)

function initsys()
    #sites=Matrix{Site}(undef,Ldim,Ldim)
    for i=1:Ldim*Ldim
        if mod(i, Ldim)==1
            v1=i+Ldim-1     #v1=-1 for Non-PBC
        else
            v1=i-1
        end
        if i>=Ldim*(Ldim-1)+1
            v2= mod(i, Ldim)==0 ? Ldim : mod(i, Ldim);
        else
            v2=i+Ldim
        end
        if mod(i, Ldim)==0
            v3=i-Ldim+1     #v1=-1 for Non-PBC
        else
            v3=i+1
        end
        if i<=Ldim
            v4= Ldim*(Ldim-1)+i;
        else
            v4=i-Ldim
        end

        global sites[i]=Site(i,rand([-1,1]),[v1 v2 v3 v4],0)
    end
end

function hamil()
    en=0;
    for i=1:Ldim*Ldim, j=1:4
        en=en+J*sites[i].spin*sites[sites[i].vec[j]].spin
    end
    return en/2
end

function dH(i, spinOLD, spinNEW)
    en=0;
    for j=1:4
        en=en+J*(spinNEW-spinOLD)*sites[sites[i].vec[j]].spin;
    end
    return en;
end

function MC(T)
    s=rand(collect(1:Ldim*Ldim));
    sspin=sites[s].spin;
    if sspin>0
        sites[s].spin=-1
    else
        sites[s].spin=1
    end
    cluster=[s];
    stack=[s];
    for i=1:4
        if sites[sites[s].vec[i]].spin==sspin
            r=rand()
            if r<1-exp(2*J/T)
                append!(cluster, sites[s].vec[i])
                append!(stack, sites[s].vec[i])
                if sspin>0
                    sites[sites[s].vec[i]].spin=-1
                else
                    sites[sites[s].vec[i]].spin=1
                end
            end
        end
    end
    #println("stack ", stack)
    #println("cluster ", cluster)

    popfirst!(stack)
    while length(stack)>0
        s2=stack[1]
        newvec=setdiff(sites[s2].vec, cluster)
        for i=1:length(newvec)
            if sites[newvec[i]].spin==sspin
                r=rand()
                if r<1-exp(2*J/T)
                    append!(cluster, newvec[i])
                    append!(stack, newvec[i])
                    if sspin>0
                        sites[newvec[i]].spin=-1
                    else
                        sites[newvec[i]].spin=1
                    end
                end
            end
        end
        #println("stack ", stack)
        #println("cluster ", cluster)
        popfirst!(stack)
    end
end

meanspin=[0.0];
spinM=[0.0];
spinMM=[0.0];

function mSpin()
    mS=0.0
    for i=1:Ldim*Ldim
        mS=mS+sites[i].spin;
    end
    mS=mS/(Ldim*Ldim);

    append!(meanspin, mS);

end

function mnsp()

    mn=sum(meanspin[i] for i=2:snap+1)/snap
    append!(spinM, mn)
    mnn=sqrt(sum((meanspin[i]-mn)^2 for i=2:snap+1)/snap)
    append!(spinMM, mnn)
    global meanspin=[0.0]
end

function mSpin2()

    for i=1:Ldim*Ldim
        append!(xArr,sites[i].spin)
    end
end

function prog()
    T=1.6;
    counter=0
    MCsteps=snap;   #10000000;
    for i=1:MCsteps
        MC(T)
        #if i>MCsteps*0.9 && counter<snap && mod(i,10^0)==0
            mSpin2()
            mSpin()
        #    counter=counter+1
        #end
    end
    mnsp();
    while T<2.9
        T=T+0.1
        counter=0
        println(T)
        for j=1:MCsteps
            MC(T)
            #if j>MCsteps*0.9 && counter<snap && mod(j,10^0)==0
                mSpin2()
                mSpin()
            #    counter=counter+1
            #end
        end
        mnsp()
    end
end

function prog2()
    T=1.6;          #AUMENTANDO EL RANGO DE T AYUDO!!!!
    counter=0
    MCsteps=10^3;   #10000000;
    for i=1:MCsteps
        MC(T)
        if i>MCsteps*0.9 && counter<snap && mod(i,10^0)==0
            mSpin2()
            mSpin()
            counter=counter+1
        end
    end
    mnsp();
    while T<2.9
        T=T+0.1
        counter=0
        println(T)
        for j=1:MCsteps
            MC(T)
            if j>MCsteps*0.9 && counter<snap && mod(j,10^0)==0
                mSpin2()
                mSpin()
                counter=counter+1
            end
        end
        mnsp()
    end
end

function main2()
    initsys()
    global meanspin=[0.0];
    global spinM=[0.0];
    global spinMM=[0.0];
    prog2();
    popfirst!(spinM)
    popfirst!(spinMM)

end

xArr=[0.0];
main2()
popfirst!(xArr)

plot(scatter(x=1.6:0.1:2.9 ,y=spinM,mode="markers"))

TRange=length(spinM)

xRes=reshape(xArr,Int(Ldim*Ldim),Int(length(xArr)/(Ldim*Ldim)))'
mean2=[xRes[j,i]-sum(xRes[:,i])/(TRange*snap) for j=1:(TRange*snap),i=1:Ldim*Ldim]
#mean=[xRes[j,i]-sum(xRes[j,:])/(Ldim*Ldim) for j=1:(TRange*snap),i=1:Ldim*Ldim]


#m=mean*transpose(mean)/((Ldim*Ldim-1))

m=transpose(mean2)*mean2/((TRange*snap-1))

d,v = eigs(m, nev=15)
#size(m)[1]

dNorm=d/sum(d)

pltEigen=plot(scatter(x=collect(1:length(dNorm)), y=dNorm, mode="markers"))
pltEigen=plot(dNorm, seriestype=:scatter,xaxis = ("ranking", :linear,  font(15, "Courier")),yaxis = ("eigenvalues", (0.0001,1),  :log,  font(15, "Courier")))
#vv=transpose(mean)*v
vv=mean2*v
plot(scatter(x=collect(1:length(v[:,1])),y=abs.(v[:,1])))

y1=vv[:,1]
y2=vv[:,2]
#plot(y1,y2, seriestype=:scatter, markercolor=append!(append!([RGB(0,0.0,1) for i=1:800],[RGB(0,1,0) for i=801:900]), [RGB(1,0.0,0) for i=901:1600]))

plot(scatter(x=y1,y=y2,mode=:markers, marker_color=[RGB(r,0,1-r) for r=0:1/(TRange*snap-1):1]))

layout=Layout(;title="Magnetization",
                     xaxis=attr(title="Temperature", showgrid=false, zeroline=false),
                     yaxis=attr(title="m", zeroline=false))
plot(scatter(x=collect(1.6:(2.9-1.6)/1400:2.9),y=abs.(y1/40), mode=:markers, marker_size="10", marker_color="rgb(142, 124, 195)"), layout)
#plot(range(1.6, stop=2.9,length=1400),abs.(y1/40))

plt=plot(y1,y2,seriestype=:scatter, markercolor=append!(append!([RGB(1,0,0) for r=1:6*(size(xRes)[1])/14],[RGB(0,1,0) for r=6*(size(xRes)[1])/14+1:8*(size(xRes)[1])/14]),[RGB(0,0,1) for r=8*(size(xRes)[1])/14+1:size(xRes)[1]]))

savefig(plt,"PCAIsing.png")
savefig(pltEigen,"EigenIsing.png")

writedlm("IsingData.dat",xRes)
#####################################################

F=svd(mean)

plot(F.S/sum(F.S), seriestype=:scatter,xaxis = ("my label", (0,16), 0:16, :linear,  font(15, "Courier")),yaxis = ("lambda", (0.0001,1),  :log,  font(15, "Courier")))

newVp=Diagonal(vcat([1,1],[0 for i=1:TRange*snap-2]))*F.Vt

Yy=transpose(mean)*newVp

plot(Yy[1,:],Yy[2,:],seriestype=:scatter)

writedlm("./Ising/x.dat",xRes)
writedlm("./Ising/M.dat",m)

xRes=readdlm("./Ising/x.dat")
