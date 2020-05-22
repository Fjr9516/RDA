function Y = Construct_Y(gnd,num_l,nClass)
Y = zeros(nClass,length(gnd));
for i = 1:num_l
    for j = 1:nClass
        if j == gnd(i)
            Y(j,i) = 1;
        else
            Y(j,i) = -1;
        end  
    end
end
end