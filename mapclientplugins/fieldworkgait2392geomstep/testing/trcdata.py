"""
Classes and functions for reading TRC motion capture files.
TRCData class by Hugh Sorby, taken from
https://github.com/mapclient-plugins/trcsourcestep
"""

class TRCData(dict):

    def load(self, filename):

        with open(filename) as f:
            contents = f.readlines()
            line_count = 0
            for line in contents:
                line_count += 1
                line = line.strip()
                if line_count == 1:
                    # File Header 1
                    sections = line.split('\t')
                    self[sections[0]] = sections[1]
                    self['DataFormat'] = sections[2]
                    data_format_count = len(sections[2].split('/'))
                    self['FileName'] = sections[3]
                elif line_count == 2:
                    # File Header 2
                    file_header_keys = line.split('\t')
                elif line_count == 3:
                    # File Header 3
                    file_header_data = line.split('\t')
                    if len(file_header_keys) == len(file_header_data):
                        for index, key in enumerate(file_header_keys):
                            if key == 'Units':
                                self[key] = file_header_data[index]
                            else:
                                self[key] = float(file_header_data[index])
                    else:
                        raise IOError('File format invalid: File header keys count (%d) is not equal to file header data count (%d)' % (len(file_header_keys), len(file_header_data)))
                elif line_count == 4:
                    # Data Header 1
                    data_header_labels = line.split('\t')
                    if data_header_labels[0] != 'Frame#':
                        raise IOError('File format not valid data header does not start with "Frame#".')
                    if data_header_labels[1] != 'Time':
                        raise IOError('File format not valid data header in position 2 is not "Time".')

                    self['Frame#'] = []
                    self['Time'] = []
                elif line_count == 5:
                    # Data Header 1
                    data_header_sublabels = line.split('\t')
                    if len(data_header_labels) != len(data_header_sublabels):
                        raise IOError('File format invalid: Data header labels count (%d) is not equal to data header sub-labels count (%d)' % (len(data_header_labels), len(data_header_sublabels)))

                    labels = []
                    for label in data_header_labels:
                        label = label.strip()
                        if len(label):
                            self[label] = []
                            labels.append(label)

                    self['Labels'] = labels
                elif line_count == 6:
                    # Blank line
                    header_read_successfully = True
                else:
                    # Data section
                    if header_read_successfully:
                        sections = line.split('\t')

                        frame = int(sections.pop(0))
                        self['Frame#'].append(frame)

                        time = float(sections.pop(0))
                        self['Time'].append(time)

                        len_section = len(sections)
                        if len_section % data_format_count == 0:
                            data = [[float(sections[subindex]) for subindex in xrange(index, index + data_format_count)] for index in xrange(0, len_section, data_format_count)]
                            self[frame] = (time, data)

                            for index, label_data in enumerate(data):
                                # Add two to the index as we want to skip 'Frame#' and 'Time'
                                self[labels[index + 2]] += [label_data]

                        else:
                            raise IOError('File format invalid: data frame %d does not match the data format' % len_section)

    def get_frame(self, frame):
        """
        Returns a dictionary of landmarks and coordinates at a particular
        frame.
        """

        landmarksNames = self['Labels']
        time, landmarksCoords = self[frame]
        landmarksNamesData = [frame, time] + landmarksCoords
        landmarks = dict(zip(landmarksNames, landmarksNamesData))
        if 'Frame#' in landmarks:
            del landmarks['Frame#']
        if 'Time' in landmarks:
            del landmarks['Time']

        return landmarks