function F = Gauss3D(x,data)

    F = x(1).*...
        exp(-(((data(:,:,:,1) - x(5)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(6)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(7)) .^2) ./ (2 * x(3) .^2))) + x(4);


end
