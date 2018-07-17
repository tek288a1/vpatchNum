function Z = cplxload(fid)
    Z = load(fid);
    Z = complex( Z(:,1), Z(:,2) );
end

