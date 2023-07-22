%regular Gallagher parity check
%konbaar@hotmail.com
%thank you for your interest

% Condition of H
% 1. (n/row)>col
% 2. col>=3
% 3. row>col
% 4. row mod n must be zero
function H = genH_regularGallagher(n, row, col)
%initial condition
% n=24
% row=6;
% col=3;

%number of message bit
k=(n/row)*col
group=n/row;

H=zeros(k,n);


%random 1'st for fisrt group

k1=0;
for m2=1:col
    
    for m1=1:row
        %
        a=0;
        while a<1,
        temp=(n/row)*rand;
        temp=(n/row)*rand;
        temp=(n/row)*rand;
        a=temp;
        end 
    
        a=int64(a+k1);
        H(a,m1)=1;
    
    end
    k1=k1+group;
end


m3=(n/row)-1;
for m1=(row+1):n
    %
    for m2=1:k
        %
        if H(m2,m1-row)==1
            %
            if mod(m2,(n/row))==0
                %
                H(m2-m3,m1)=1;
            else
                %
                H(m2+1,m1)=1;
            end           
        end
    end
end 




%follow Gallagher's pattern
for m1=1:(group)
    %
    H(m1,:)=0;
end

k1=1;
k2=0;
for m2=1:(group)

    for m1=1:row
        %
        H(k1,m1+k2)=1;
    end

    k1=k1+1;
    k2=k2+row;

end

%check for constant column weight
Wcol=0;
for c=1:n
    %
    for r=1:k
        %
        if H(r,c)==1
            Wcol=Wcol+1;
        end
    end
end

if Wcol==(col*n)
    %
    fprintf('\nOk');
    column_weight=col
else
    %
    fprintf('\n I am not ok');
end


%check for constant row weight
Wrow=0;
for r=1:k
    %
    for c=1:n
        %
        if H(r,c)==1
            Wrow=Wrow+1;
        end
    end
end

if Wrow==(row*k)
    %
    fprintf('\nOk');
    row_weight=row
else
    %
    fprintf('\n I am not ok');
end

%Display Parity Check
H
LDPC_girth4a(H)