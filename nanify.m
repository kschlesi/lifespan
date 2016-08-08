function xmat = nanify(xmat,dr)

xmat(xmat==dr) = NaN;

end