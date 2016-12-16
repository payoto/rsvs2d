xDoc = xmlread('C:\Users\am14795\Chrome Local Downloads\takeout-20161216T142451Z\Takeout\Location History\LocationHistory.kml')
gxcoord = xDoc.getElementsByTagName('gx:coord');
when = xDoc.getElementsByTagName('when');

len = gxcoord.getLength;
for i = 1:len
    disp(i)
    elem = gxcoord.item(i);
    bf = elem.getFirstChild.getData;
    bf = sscanf(char(bf),'%f %f %f');
    data(i,1:2) = bf(1:2);
end