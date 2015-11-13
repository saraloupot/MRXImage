function dist = gaussianDistance(source,sourcelocation,sourcestrength)
load('headerVariables');
xsolution = interplocation(sourcelocation,sourcestrength,points);
xsolution = xsolution./(sum(xsolution));
solution = reshape(xsolution,ndimx,ndimy,ndimz);


reconlocation = source(:,1:3);
reconstrength = zeros(size(reconlocation));
reconstrength(:,3) = source(:,4);
xrecon = interplocation(reconlocation,reconstrength, points);
xrecon = xrecon./(sum(xrecon));
recon = reshape(xrecon,ndimx,ndimy,ndimz);

blurredReconImage = smooth3(recon,'gaussian',[2*ndimx+1,2*ndimy+1,2*ndimz+1],2);
blurredSourceImage = smooth3(solution,'gaussian',[2*ndimx+1,2*ndimy+1,2*ndimz+1],2);
dist=norm(blurredReconImage(:)-blurredSourceImage(:));

function v = interplocation(sourcelocation,sourcestrength,points)
vi = [];
for ii = 1:size(sourcelocation,1)
    [xi yi zi]=ndgrid(-.5:.5:.5,-.5:.5:.5,-.5:.5:.5);
    xi=xi+sourcelocation(ii,1);
    yi = yi+sourcelocation(ii,2);
    zi = zi+sourcelocation(ii,3);
    xs = zeros(size(xi));
    xs(2,2,2)=sourcestrength(ii,3);
    vi(:,ii) = griddata(xi(:),yi(:),zi(:),xs(:),points(:,1),points(:,2),points(:,3),'nearest');
end
v=sum(vi,2);