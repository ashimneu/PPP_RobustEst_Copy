function error_val = err(p,xk)
error_val = norm(p.Grdpos.pos(:,1)-xk(1:3));
end

