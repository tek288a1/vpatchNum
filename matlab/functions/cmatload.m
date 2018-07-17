function M = cmatload(fid)
    Z = load(fid);
    Z = complex( Z(:,1), Z(:,2) );
    n = sqrt(length(Z));
    M = reshape(Z, n, n);
end

