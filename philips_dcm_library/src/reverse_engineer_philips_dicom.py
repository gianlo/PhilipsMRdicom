import gzip
import dicom


of = gzip.GzipFile(r'D:\var\common\T2.phoetal.dynamic.sag.dcm.metadata.gz')
lines = of.readlines()
of.close()

frametimings = dict()
for l in lines:
	if l.startswith('#'):
		continue
	if not l.find(' = ')>-1:
		continue
	tag, value = tuple(l.split(' = '))
	
	if tag.find('(2005,100b)')>-1:		
		frame_number = tag[len('PerFrameFunctionalGroupsSequence['):tag.find(']')]
		key = int(frame_number)-1
		if key not in frametimings:
			frametimings[key] = [0, 0.]
		frametimings[key][1] = float(value)
	elif tag.find('].InStackPositionNumber')>-1:
		frame_number = tag[len('PerFrameFunctionalGroupsSequence['):tag.find(']')]
		key = int(frame_number)-1
		if key not in frametimings:
			frametimings[key] = [0, 0.]
		frametimings[key][0] = int(value)

ft = zip(map(int,frametimings.keys()), frametimings.values())
z = sorted(ft, key=lambda x:(x[1][1],x[1][0]))
df = dicom.read_file(r'D:\var\common\T2.phoetal.dynamic.sag.dcm')
