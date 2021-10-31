#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Build the Apollo data archive for the PDS

Note that this requires mseed2ascii code. See note below. 
'''


from lxml import etree
from obspy import read
import os
import fnmatch
import shutil
import math
import pandas as pd
import re
from obspy.core.utcdatetime import UTCDateTime
from obspy import read_inventory
from subprocess import Popen, PIPE
from datetime import datetime
import logging
from copy import deepcopy
from collections import OrderedDict 

# when using find, can use wildcards like this 
# see https://lxml.de/2.2/tutorial.html
# - all tags in a namespace: "{namespace}*"
# - a local name 'tag' in any (or no) namespace: "{*}tag"
# - a tag without namespace: "{}tag"
# - all tags without namespace: "{}*"


logging.basicConfig(filename='/Users/cnunn/Downloads/log.txt', filemode='w', level=logging.DEBUG)

LOGICAL_IDENTIFIER_ROOT = 'urn:nasa:pds:apollo_pse:data_seed'
MISSION_PHASE_NAME = 'SURFACE MISSION'
LID_REFERENCE_ROOT = 'urn:nasa:pds:apollo_pse:data_seed'
LIDVID_REFERENCE_ROOT = 'urn:nasa:pds:apollo_pse:data_table'
LID_REFERENCE_CONTEXT_ROOT = 'urn:nasa:pds:context'


LOGICAL_IDENTIFIER_ROOT_GEOCSV = 'urn:nasa:pds:apollo_pse:data_table'
# MISSION_PHASE_NAME = 'SURFACE MISSION'
LID_REFERENCE_ROOT_GEOCSV = 'urn:nasa:pds:apollo_pse:data_table:stationxml'
LIDVID_REFERENCE_ROOT_GEOCSV = 'urn:nasa:pds:apollo_pse:data_seed'

TEST=False # set to True to test a single file 

DELTA = 0.1509433962


# c ------------ for Release 1 ----------------
# c XA dataless.xa.2021.060.1.seed

DATALESS_DICTIONARY = {
     '1' : {
    'S11': 'xa.0',
    'S12': 'xa.0',
    'S14': 'xa.0',
    'S15': 'xa.0',
    'S16': 'xa.0',
    }
}

def load_template(name):
    ''' Takes an xml file as input. Outputs ElementTree and element'''
    # specify parser setting
    parser = etree.XMLParser(strip_cdata=False,remove_blank_text=True)
    # pass parser to do the actual parsing
    tree = etree.parse(name, parser)

    root = tree.getroot()
    return tree, root

def make_seed_xml(output_dir,seed_basename,release_number,seed_col):
    # Make the seed label

    # read in the miniseed template
    template_tree, product_native = load_template('PDS_files/miniseed_template.xml')

    # xa.s12..mhz.1976.061.0.mseed
    net_upper, sta_upper, loc_upper, chan_upper, year, jday, rev, _ = seed_basename.upper().split('.')

    net, sta, loc, chan, year, jday, rev, _ = seed_basename.lower().split('.')

    sta_number = sta.lower()[1:3]

    seed_filename = os.path.join(output_dir,seed_basename)

    file_modification_date = UTCDateTime(os.path.getmtime(seed_filename))

    # read in the seed file 
    st = read(seed_filename)
    st.merge()
    for tr in st:     
        npts = tr.stats.npts
        starttime = tr.stats.starttime
        endtime = tr.stats.endtime
        sampling_rate_Hz = tr.stats.sampling_rate
        break

    # test one, replace with new version
    # net = 'xb'
    # sta = 'elyhk'
    # loc = '45'
    # chan = 'uea'
    # year = '2019'
    # jday = '273'
    # rev = '4'

    # TODO - fix 

    version_id = '1.0'
    seed_col = update_collection_seed_file(seed_col,net, sta, loc, chan, year, jday, rev, version_id)

    # get the namespace for apollo
    namespace = template_tree.xpath('namespace-uri(.)')


    # loop through the template tree
    for identification_area in product_native.findall('{%s}Identification_Area' % namespace):

        logical_identifier = identification_area.find('{%s}logical_identifier' % namespace)
        # urn:nasa:pds:apollo_pse:data_seed:<!-- |net.sta.loc.chan.year.jday.rev| -->

        logical_identifier.text = '{}:{}.{}.{}.{}.{}.{}.{}'.format(
          LOGICAL_IDENTIFIER_ROOT, net, sta, loc, chan, year, jday, rev)

        # identifier_text is required for the collection files
        identifier_text = 'P,{}::{}'.format(logical_identifier.text,version_id)

        ver_id = identification_area.find('{%s}version_id' % namespace)
        ver_id.text = version_id

        for modication_history in identification_area.findall('{%s}Modification_History' % namespace):
            for modification_detail in modication_history.findall('{%s}Modification_Detail' % namespace):
                modification_date = modification_detail.find('{%s}modification_date' % namespace)
                print(modification_date.text)
                # <!-- |@CurrentUTCDateTime| -->Z
                modification_date.text = file_modification_date.strftime("%Y-%m-%dZ")


                ver_id = modification_detail.find('{%s}version_id' % namespace)
                ver_id.text = version_id

    for context_area in product_native.findall('{%s}Context_Area' % namespace):
        for time_coordinates in context_area.findall('{%s}Time_Coordinates' % namespace):
            start_date_time=time_coordinates.find('{%s}start_date_time' % namespace)
            # <!-- |start_date_time| -->Z
            start_date_time.text = str(starttime)
            stop_date_time=time_coordinates.find('{%s}stop_date_time' % namespace)
            # <!-- |stop_date_time| -->Z
            stop_date_time.text = str(endtime)

        investigation_area = context_area.find('{*}Investigation_Area')
        name = investigation_area.find('{*}name')
        name.text = 'APOLLO {}'.format(sta_number)

        internal_reference = investigation_area.find('{*}Internal_Reference')
        lid_reference = internal_reference.find('{*}lid_reference')
        lid_reference.text = '{}:investigation:mission.apollo_{}'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)

        observing_system = context_area.find('{*}Observing_System')
        for observing_system_component in observing_system.findall('{%s}Observing_System_Component' % namespace):
            name = observing_system_component.find('{*}name')
            type1 = observing_system_component.find('{*}type')

            internal_reference = observing_system_component.find('{*}Internal_Reference')
            lid_reference = internal_reference.find('{*}lid_reference')

            if type1.text == 'Host':
                name.text = 'APOLLO {} LUNAR SURFACE EXPERIMENTS PACKAGE'.format(sta_number)
                if sta_number == '11':
                    lid_reference.text = '{}:instrument_host:spacecraft.a{}e'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
                else:
                    lid_reference.text = '{}:instrument_host:spacecraft.a{}a'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
            else: 
                name.text = 'Apollo {} Passive Seismic Experiment (PSE)'.format(sta_number)
                if sta_number == '11':
                    lid_reference.text = '{}:instrument:pse.a{}e'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
                else:
                    lid_reference.text = '{}:instrument:pse.a{}a'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
        # YYYY
                                                                 
        mission_area = context_area.find('{*}Mission_Area')

        seismic_parameters = mission_area.find('{*}Seismic_Parameters')
        metadata_location = seismic_parameters.find('{*}Metadata_Location')
        metadata_file_name = metadata_location.find('{*}metadata_file_name')

        # <apollo:metadata_file_name>dataless.xa.2021.060.seed</apollo:metadata_file_name>
        metadata_file_name.text = 'dataless.{}.seed'.format(DATALESS_DICTIONARY[release_number][sta.upper()])

        meta_internal_reference = metadata_location.find('{*}Internal_Reference')
        lid_reference = meta_internal_reference.find('{*}lid_reference')
        # <lid_reference>urn:nasa:pds:apollo_pse:data_seed:dataless.xa.2021.001</lid_reference>
        lid_reference.text = '{}:dataless.{}'.format(LID_REFERENCE_ROOT, DATALESS_DICTIONARY[release_number][sta.upper()])

        ASCII_Equivalent = seismic_parameters.find('{*}ASCII_Equivalent')
        ascii_equivalent_file_name = ASCII_Equivalent.find('{*}ascii_equivalent_file_name')
        # <apollo:ascii_equivalent_file_name>xa.s16.01.afr.1976.066.1.1.a.csv</apollo:ascii_equivalent_file_name>
        ascii_equivalent_file_name.text = '{}.{}.{}.{}.{}.{}.{}.a.csv'.format(net,sta,loc,chan,year,jday,rev)

        internal_reference = ASCII_Equivalent.find('{*}Internal_Reference')
        lidvid_reference = internal_reference.find('{*}lidvid_reference')

        # <lidvid_reference>urn:nasa:pds:apollo_pse:data_table:xa.s16.01.afr.1976.066.1.1.a::1.0</lidvid_reference>
        lidvid_reference.text = '{}:{}.{}.{}.{}.{}.{}.{}.a::{}'.format(
                      LIDVID_REFERENCE_ROOT,net,sta,loc,chan,year,jday,rev,version_id)

        station_tag = seismic_parameters.find('{*}station')
        station_tag.text = sta_upper
        
        channel_tag = seismic_parameters.find('{*}channel')
        channel_tag.text = chan_upper

        location_tag = seismic_parameters.find('{*}location')
        if loc == '':
            seismic_parameters.remove(location_tag)
        else:
            location_tag.text = loc
        
        sample_count = seismic_parameters.find('{*}sample_count')
        sample_count.text = str(npts)

        sampling_rate = seismic_parameters.find('{*}sampling_rate')
        sampling_rate.text = str(sampling_rate_Hz)

    for file_area_native in product_native.findall('{%s}File_Area_Native' % namespace):
        for file1 in file_area_native.findall('{%s}File' % namespace):
            file_name=file1.find('{%s}file_name' % namespace)
            # <!-- |net.sta.loc.chan.year.jday.rev| -->.mseed
            file_name.text = '{}.{}.{}.{}.{}.{}.{}.mseed'.format(net,sta,loc,
              chan,year,jday,rev)

            creation_date_time=file1.find('{%s}creation_date_time' % namespace)
            # <!-- |creation_date_time| -->
            creation_date_time.text = str(file_modification_date)

    # create a new XML file with the results
    xml_data = etree.tostring(template_tree, xml_declaration=True, encoding="utf-8", pretty_print=True).decode("utf-8")

#     # product native line not being formatted nicely, so replace it 
    repl = '''<Product_Native
      xmlns="http://pds.nasa.gov/pds4/pds/v1"
      xmlns:apollo="http://pds.nasa.gov/pds4/mission/apollo/v1"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://pds.nasa.gov/pds4/pds/v1
      https://pds.nasa.gov/pds4/pds/v1/PDS4_PDS_1G00.xsd
      http://pds.nasa.gov/pds4/mission/apollo/v1
      https://pds.nasa.gov/pds4/mission/apollo/v1/PDS4_APOLLO_1G00_1000.xsd">
'''
    pattern = '<Product_Native.*'
    xml_data = re.sub(pattern, repl, xml_data, count=0, flags=0)

    base = os.path.splitext(seed_filename)[0]
    xml_filename = base + '.xml'

    # write the seed label xml file 
    xml_file = open(xml_filename, 'w', newline='\r\n')
    print('File to write ', xml_filename)
    xml_file.write(xml_data)
    xml_file.close()

def make_dataless_xml(output_dir,dataless_basename,release_number,dataless_col):
    # Make the seed label

    # read in the miniseed template
    template_tree, product_native = load_template('PDS_files/dataless_template.xml')

    # dataless.xa.2021.060.1.seed
    # dataless.xa.2021.060.1.xml
    # filetype, net_upper, sta_upper, loc_upper, chan_upper, year, jday, rev, _ = dataless_basename.upper().split('.')
    print(dataless_basename)
#     exit()
# dataless.xa.1969-1977.0.dataless

    _, net, rev, _ = dataless_basename.lower().split('.')

    dataless_filename = os.path.join(output_dir,dataless_basename)

    file_modification_date = UTCDateTime(os.path.getmtime(dataless_filename))

    # test one, replace with new version
    # net = 'xb'
    # sta = 'elyhk'
    # loc = '45'
    # chan = 'uea'
    # year = '2019'
    # jday = '273'
    # rev = '4'


    version_id = '1.0'


    dataless_col = update_collection_dataless_file(dataless_col,net, rev, version_id)

    # get the namespace for apollo
    namespace = template_tree.xpath('namespace-uri(.)')

    # loop through the template tree
    for identification_area in product_native.findall('{%s}Identification_Area' % namespace):
    
        logical_identifier = identification_area.find('{%s}logical_identifier' % namespace)
    
        # identifier_text is required for the collection files
        identifier_text = 'P,{}::{}'.format(logical_identifier.text,version_id)
    
        ver_id = identification_area.find('{%s}version_id' % namespace)
        ver_id.text = version_id
    
        for modication_history in identification_area.findall('{%s}Modification_History' % namespace):
            for modification_detail in modication_history.findall('{%s}Modification_Detail' % namespace):
                modification_date = modification_detail.find('{%s}modification_date' % namespace)
                # <!-- |@CurrentUTCDateTime| -->Z
                modification_date.text = file_modification_date.strftime("%Y-%m-%dZ")
    
                ver_id = modification_detail.find('{%s}version_id' % namespace)
                ver_id.text = version_id

            context_area = product_native.find('{*}Context_Area')
            mission_area = context_area.find('{*}Mission_Area')
            seismic_parameters = mission_area.find('{*}Seismic_Parameters')
            ASCII_Equivalent = seismic_parameters.find('{*}ASCII_Equivalent')

            internal_reference = ASCII_Equivalent.find('{*}Internal_Reference')
            lidvid_reference = internal_reference.find('{*}lidvid_reference')


    for file_area_native in product_native.findall('{%s}File_Area_Native' % namespace):
        for file1 in file_area_native.findall('{%s}File' % namespace):
            creation_date_time=file1.find('{%s}creation_date_time' % namespace)
            creation_date_time.text = str(file_modification_date)

    # create a new XML file with the results
    xml_data = etree.tostring(template_tree, xml_declaration=True, encoding="utf-8", pretty_print=True).decode("utf-8")



    # Product_Native line not being formatted nicely, so replace it 
    repl = '''<Product_Native
      xmlns="http://pds.nasa.gov/pds4/pds/v1"
      xmlns:apollo="http://pds.nasa.gov/pds4/mission/apollo/v1"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://pds.nasa.gov/pds4/pds/v1
      https://pds.nasa.gov/pds4/pds/v1/PDS4_PDS_1G00.xsd
      http://pds.nasa.gov/pds4/mission/apollo/v1
      https://pds.nasa.gov/pds4/mission/apollo/v1/PDS4_APOLLO_1G00_1000.xsd">
'''

    pattern = '<Product_Native.*'
    xml_data = re.sub(pattern, repl, xml_data, count=0, flags=0)

    base = os.path.splitext(dataless_filename)[0]
    xml_filename = base + '.xml'

    # write the dataless label xml file 
    xml_file = open(xml_filename, 'w', newline='\r\n')
    print('File to write ', xml_filename)

    xml_file.write(xml_data)



def make_stationxml_xml(output_dir,stationxml_basename,release_number,stationxml_col):
    # Make the seed label

    # read in the miniseed template
    template_tree, product_native = load_template('PDS_files/stationxml_template.xml')


    print(stationxml_basename)
    # stationxml.xa.1969-1977.0.sxml
    filetype, net, rev, _ = stationxml_basename.lower().split('.')

    stationxml_filename = os.path.join(output_dir,stationxml_basename)

    file_modification_date = UTCDateTime(os.path.getmtime(stationxml_filename))


    # test one, replace with new version
    # net = 'xb'
    # sta = 'elyhk'
    # loc = '45'
    # chan = 'uea'
    # year = '2019'
    # jday = '273'
    # rev = '4'

    # P,urn:nasa:pds:apollo_pse:data_seed:dataless.xb.elyse.2020.062::1.0
    version_id = '1.0'
    stationxml_col = update_collection_stationxml_file(stationxml_col,net,rev,version_id)

    # get the namespace for insight
    namespace = template_tree.xpath('namespace-uri(.)')

    # loop through the template tree
    for identification_area in product_native.findall('{%s}Identification_Area' % namespace):
    
        logical_identifier = identification_area.find('{%s}logical_identifier' % namespace)
        # urn:nasa:pds:apollo_pse:data_seed:<!-- |net.sta.loc.chan.year.jday.rev| -->
    
        # logical_identifier.text = '{}:stationxml.{}.{}'.format(
        #   LIDVID_REFERENCE_ROOT, net, rev)
    
        # identifier_text is required for the collection files
        identifier_text = 'P,{}::{}'.format(logical_identifier.text,version_id)
    
        ver_id = identification_area.find('{%s}version_id' % namespace)
        ver_id.text = version_id
    
        for modication_history in identification_area.findall('{%s}Modification_History' % namespace):
            for modification_detail in modication_history.findall('{%s}Modification_Detail' % namespace):
                modification_date = modification_detail.find('{%s}modification_date' % namespace)
                # <!-- |@CurrentUTCDateTime| -->Z
                modification_date.text = file_modification_date.strftime("%Y-%m-%dZ")
    
                ver_id = modification_detail.find('{%s}version_id' % namespace)
                ver_id.text = version_id

    for context_area in product_native.findall('{%s}Context_Area' % namespace):
        for mission_area in context_area.findall('{%s}Mission_Area' % namespace):
            seismic_parameters = mission_area.find('{*}Seismic_Parameters')
            SEED_equivalent = seismic_parameters.find('{*}SEED_Equivalent')
            # seed_file_name = SEED_equivalent.find('{*}seed_file_name')
            # <apollo:seed_file_name>dataless.xa.2021.001.seed</apollo:seed_file_name>
            # seed_file_name.text = 'dataless.{}.{}.seed'.format(net,rev)

            internal_reference = SEED_equivalent.find('{*}Internal_Reference')
            lidvid_reference = internal_reference.find('{*}lidvid_reference')
            # <lidvid_reference>urn:nasa:pds:apollo_pse:data_seed:dataless.xa.99::1.0</lidvid_reference>
            # lidvid_reference.text = '{}:dataless.{}.{}::{}'.format(
            #   LOGICAL_IDENTIFIER_ROOT,net,rev,version_id)

    for file_area_ancillary in product_native.findall('{%s}File_Area_Ancillary' % namespace):
        
        for file1 in file_area_ancillary.findall('{%s}File' % namespace):
            file_name=file1.find('{%s}file_name' % namespace)
            # dataless.xa.2021.060.1.seed
            # file_name.text = 'stationxml.{}.{}.sxml'.format(net,rev)
    # 
            creation_date_time=file1.find('{%s}creation_date_time' % namespace)
            # <!-- |creation_date_time| -->
            creation_date_time.text = str(file_modification_date)

    # create a new XML file with the results
    xml_data = etree.tostring(template_tree, xml_declaration=True, encoding="utf-8", pretty_print=True).decode("utf-8")

    # print(xml_data)
    # exit()

    # product native line not being formatted nicely, so replace it 
    repl = '''<Product_Ancillary
      xmlns="http://pds.nasa.gov/pds4/pds/v1"
      xmlns:apollo="http://pds.nasa.gov/pds4/mission/apollo/v1"
      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
      xsi:schemaLocation="http://pds.nasa.gov/pds4/pds/v1
      https://pds.nasa.gov/pds4/pds/v1/PDS4_PDS_1G00.xsd
      http://pds.nasa.gov/pds4/mission/apollo/v1
      https://pds.nasa.gov/pds4/mission/apollo/v1/PDS4_APOLLO_1G00_1000.xsd">
'''

    # XXXX
    pattern = '<Product_Ancillary.*'
    xml_data = re.sub(pattern, repl, xml_data, count=0, flags=0)

    base = os.path.splitext(stationxml_filename)[0]

    xml_filename = base + '.xml'

    # write the seed label xml file 
    xml_file = open(xml_filename, 'w', newline='\r\n')
    print('File to write ', xml_filename)

    xml_file.write(xml_data)

def make_geocsv_xml(output_dir,csv_basename,release_number,geocsv_col):

    # make the geocsv file 
    # print(csv_basename)
    # read in the geocsv template
    template_tree_geocsv, product_native_geocsv = load_template('PDS_files/geocsv_template.xml')

    # /Users/cnunn/lunar_data/PDART_PDS_output/xa/continuous_waveform/s12/1976/066/xa.s12.01.mhz.1976.066.1.1.a.csv

    net, sta, loc, chan, year, jday, rev, _, _ = csv_basename.lower().split('.')

    net_upper, sta_upper, loc_upper, chan_upper, _, _, _, _, _ = csv_basename.upper().split('.')

    sta_number = sta.lower()[1:3]

    csv_filename = os.path.join(output_dir,csv_basename)

    file_modification_date = UTCDateTime(os.path.getmtime(csv_filename))
    
    npts = []
    start_times = []
    
    with open(csv_filename) as fp:
        for line in fp:
            # check for an error
            # TODO remove this line if we have better error handling elsewhere
            if 'Error' in line:
                print(line)
                # this is a major error and we can't continue
                # remove the file with the error 
                csv_filename_str = str(csv_filename)
                if os.path.exists(csv_filename):
                    os.remove(csv_filename)
                return (line, csv_filename_str)
            if '# sample_count:' in line:
                npts.append(int(line.strip('\n').split(' ')[-1]))
            if '# sample_rate_hz:' in line:
                sample_rate_hz=line.strip('\n').split(' ')[-1]
            if '# start_time:' in line:
                start_times.append(line.strip('\n').split(' ')[-1])
            if '# field_type:' in line:
                field_type2 = line.strip('\n').split(' ')[-1]
        end_time = line.split(',')[0]
        npts_total = str(sum(npts))


    # # if the date is earlier than REVISIONS_PRIOR_TO_DATE, then it is a revision
    # # TODO what if there is more than one revision
    # if UTCDateTime(start_times[0]) < UTCDateTime(REVISIONS_PRIOR_TO_DATE):
    #     version_id = '2.0'
    # else:
    #     version_id = '1.0'

    version_id = '1.0'
    geocsv_col = update_collection_geocsv_file(geocsv_col,net, sta, loc, chan, year, jday, rev, version_id)

    offsets = []
    object_lengths = []
    col_header_offsets = []
    col_header_object_length = '14' # the byte length of Time, Sample and its carriage returns
    table_offsets = []

    offset_test = bytes('# dataset', 'utf-8')
    offset_test_col = bytes('Time, Sample', 'utf-8')
    # reread the file to get the byte offsets
    with open(csv_filename, 'rb') as fp:
        s = fp.read()
        index = 0
        while index < len(s):
            index = s.find(offset_test, index)
            if index == -1:
                break
            offsets.append(str(index))

            index += len(offset_test)
            index = s.find(offset_test_col, index)

            if index == -1:
                break
            col_header_offsets.append(str(index))

    # calculate byte length (object_length) from offsets - col_header_offsets
    for (a, b) in zip(col_header_offsets, offsets):
        object_lengths.append(str(int(a)-int(b)))

    # calculate byte length (object_length) from col_header_object_length 
    # plus the length of the column header
    for a in col_header_offsets:
        table_offsets.append(str(int(a)+int(col_header_object_length)))


    # print('object_lengths', object_lengths)
    # print('table_offsets', table_offsets)

    # print('temp !!!!')
    # exit()

    # get the insight namespace
    namespace = template_tree_geocsv.xpath('namespace-uri(.)')

    # loop through the template file 
    for identification_area in product_native_geocsv.findall('{%s}Identification_Area' % namespace):

        logical_identifier = identification_area.find('{%s}logical_identifier' % namespace)
        # urn:nasa:pds:apollo_pse:data_table:<!-- |net.sta.loc.chan.year.jday.rev| -->.a
        # 'WARNING Is the geo csv identifier correct?'

        logical_identifier.text = '{}:{}.{}.{}.{}.{}.{}.{}.a'.format(
          LOGICAL_IDENTIFIER_ROOT_GEOCSV, net, sta, loc, chan, year, jday, rev)
        identifier_text = 'P,{}::{}'.format(logical_identifier.text,version_id)

        ver_id = identification_area.find('{%s}version_id' % namespace)
        ver_id.text = version_id
        # 1.0

        for modication_history in identification_area.findall('{%s}Modification_History' % namespace):
            for modification_detail in modication_history.findall('{%s}Modification_Detail' % namespace):
                modification_date = modification_detail.find('{%s}modification_date' % namespace)
                # <!-- |@CurrentUTCDateTime| -->Z
                modification_date.text = file_modification_date.strftime("%Y-%m-%dZ")

                ver_id = modification_detail.find('{%s}version_id' % namespace)
                # 1.0
                ver_id.text = version_id

    for context_area in product_native_geocsv.findall('{%s}Observation_Area' % namespace):
        for time_coordinates in context_area.findall('{%s}Time_Coordinates' % namespace):
            start_date_time=time_coordinates.find('{%s}start_date_time' % namespace)
            # <!-- |start_date_time| -->Z
            start_date_time.text = start_times[0]
            stop_date_time=time_coordinates.find('{%s}stop_date_time' % namespace)
            # <!-- |stop_date_time| -->Z
            stop_date_time.text = end_time

        investigation_area = context_area.find('{*}Investigation_Area')
        name = investigation_area.find('{*}name')
        name.text = 'APOLLO {}'.format(sta_number)

        internal_reference = investigation_area.find('{*}Internal_Reference')
        lid_reference = internal_reference.find('{*}lid_reference')
        lid_reference.text = '{}:investigation:mission.apollo_{}'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)

        observing_system = context_area.find('{*}Observing_System')
        for observing_system_component in observing_system.findall('{%s}Observing_System_Component' % namespace):
            name = observing_system_component.find('{*}name')
            type1 = observing_system_component.find('{*}type')

            internal_reference = observing_system_component.find('{*}Internal_Reference')
            lid_reference = internal_reference.find('{*}lid_reference')

            if type1.text == 'Host':
                name.text = 'APOLLO {} LUNAR SURFACE EXPERIMENTS PACKAGE'.format(sta_number)
                if sta_number == '11':
                    lid_reference.text = '{}:instrument_host:spacecraft.a{}e'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
                else:
                    lid_reference.text = '{}:instrument_host:spacecraft.a{}a'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
            else: 
                name.text = 'Apollo {} Passive Seismic Experiment (PSE)'.format(sta_number)
                if sta_number == '11':
                    lid_reference.text = '{}:instrument:pse.a{}e'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
                else:
                    lid_reference.text = '{}:instrument:pse.a{}a'.format(LID_REFERENCE_CONTEXT_ROOT, sta_number)
        # YYYY
        
        for mission_area in context_area.findall('{%s}Mission_Area' % namespace):
            for elem in mission_area:
                if 'Observation_Information' in elem.tag:
                    for subelem in elem: 
                        if 'release_number' in subelem.tag: 
                            # <!-- |release_number| -->
                            subelem.text = release_number
                        if 'mission_phase_name' in subelem.tag:
                            # <!-- |mission_phase_name| -->
                            subelem.text = MISSION_PHASE_NAME

            seismic_parameters = mission_area.find('{*}Seismic_Parameters')

            metadata_location = seismic_parameters.find('{*}Metadata_Location')

            metadata_file_name = metadata_location.find('{*}metadata_file_name')
            metadata_file_name.text = 'stationxml.{}.sxml'.format(DATALESS_DICTIONARY[release_number][sta_upper])
            full_metadata_file_name = metadata_file_name.text

            internal_reference = metadata_location.find('{*}Internal_Reference')
            lid_reference = internal_reference.find('{*}lid_reference')
            lid_reference.text = '{}.{}'.format(LID_REFERENCE_ROOT_GEOCSV,DATALESS_DICTIONARY[release_number][sta_upper])


            SEED_Equivalent = seismic_parameters.find('{*}SEED_Equivalent')
            seed_file_name = SEED_Equivalent.find('{*}seed_file_name')
            seed_file_name.text = '{}.{}.{}.{}.{}.{}.{}.mseed'.format(net,sta,loc,chan,year,jday,rev)
            internal_reference = SEED_Equivalent.find('{*}Internal_Reference')
            lidvid_reference = internal_reference.find('{*}lidvid_reference')
            lidvid_reference.text = '{}:{}.{}.{}.{}.{}.{}.{}::{}'.format(
              LIDVID_REFERENCE_ROOT_GEOCSV,net,sta,loc,chan,year,jday,rev,version_id)

            station_tag = seismic_parameters.find('{*}station')
            station_tag.text = sta_upper
            
            channel_tag = seismic_parameters.find('{*}channel')
            channel_tag.text = chan_upper

            location_tag = seismic_parameters.find('{*}location')
            if loc == '':
                seismic_parameters.remove(location_tag)
            else:
                location_tag.text = loc
            
            sample_count = seismic_parameters.find('{*}sample_count')
            sample_count.text = npts_total

            sampling_rate = seismic_parameters.find('{*}sampling_rate')
            sampling_rate.text = sample_rate_hz

    for file_area_observational in product_native_geocsv.findall('{%s}File_Area_Observational' % namespace):
        for file1 in file_area_observational.findall('{%s}File' % namespace):
            file_name=file1.find('{%s}file_name' % namespace)
            # <!-- |net.sta.loc.chan.year.jday.rev| -->.a.csv
            file_name.text = '{}.{}.{}.{}.{}.{}.{}.a.csv'.format(net,sta,loc,
              chan,year,jday,rev)

            creation_date_time=file1.find('{%s}creation_date_time' % namespace)
            # <!-- |creation_date_time| -->
            creation_date_time.text = str(file_modification_date)

        # read in the two header tags, and then remove them from the template
        for header in file_area_observational.findall('{%s}Header' % namespace):
            for subelem in header:
                if subelem.text == 'Table 1 Comments':
                    header1 = deepcopy(header)
                    file_area_observational.remove(header)
                if subelem.text == 'Table 1 Column Headings':
                    header2 = deepcopy(header)
                    file_area_observational.remove(header)

        # set the column data types and descriptions, if required
        table_delimited = file_area_observational.find('{*}Table_Delimited')
        # remove table_delimited the template 
        file_area_observational.remove(table_delimited)

        record_delimited = table_delimited.find('{*}Record_Delimited')
        for field_delimited in record_delimited.findall('{%s}Field_Delimited' % namespace):  
            name = field_delimited.find('{*}name')
            if name.text == 'Sample':
                description = field_delimited.find('{*}description')
                if chan_upper == 'AGR':
                    description.text = 'Number of the ground station recording the download.'
                elif chan_upper == 'ATT':
                    description.text = 'Time the signal was recorded on Earth at the ground station. Time in seconds since the beginning of the epoch (1 January 1970). The timing can be negative for early in the mission.'
                elif chan_upper == 'AFR':
                    description.text = 'Frame number of the sample (frame numbers range from 0 to 89 and were recorded on the data sampler).'

                if field_type2 == 'FLOAT':
                    data_type = field_delimited.find('{*}data_type')
                    data_type.text = 'ASCII_Real'

        # repeat for each block 
        for i, _ in enumerate(start_times):
        
            # copy the headers and block_table_delimited
            block_header1 = deepcopy(header1)
            block_header2 = deepcopy(header2)
            block_table_delimited = deepcopy(table_delimited)
        
            # update the first header 
            for subelem in block_header1:
                if 'name' in subelem.tag:
                    subelem.text = 'Table {} Comments'.format(str(i+1))
                if 'offset' in subelem.tag:
                    subelem.text = offsets[i]
                if 'object_length' in subelem.tag:
                    subelem.text = object_lengths[i]
        
            # update the second header         
            for subelem in block_header2:
                if 'name' in subelem.tag:
                    subelem.text = 'Table {} Column Headings'.format(str(i+1))
                if 'offset' in subelem.tag:
                    subelem.text = col_header_offsets[i]
                if 'object_length' in subelem.tag:
                    subelem.text = col_header_object_length
        
            # update the block_table_delimited 
            for subelem in block_table_delimited:
                if 'name' in subelem.tag:
                    subelem.text = 'Table {} in Apollo GeoCSV File {}.{}.{}.{}.{}.{}.{}.a.csv'.format(str(i+1),net,sta,loc,
                                  chan,year,jday,rev)
                if 'offset' in subelem.tag:
                    subelem.text = table_offsets[i]
                if 'records' in subelem.tag:
                    subelem.text = records=str(npts[i])
        
            # add them to the template 
            file_area_observational.append(block_header1)
            file_area_observational.append(block_header2)
            file_area_observational.append(block_table_delimited)
            
    # create a new XML file with the results
    xml_data = etree.tostring(template_tree_geocsv, xml_declaration=True, encoding="utf-8", pretty_print=True).decode("utf-8")

    repl='''<Product_Observational
        xmlns="http://pds.nasa.gov/pds4/pds/v1"
        xmlns:apollo="http://pds.nasa.gov/pds4/mission/apollo/v1"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:schemaLocation="http://pds.nasa.gov/pds4/pds/v1
        https://pds.nasa.gov/pds4/pds/v1/PDS4_PDS_1G00.xsd
        http://pds.nasa.gov/pds4/mission/apollo/v1
        https://pds.nasa.gov/pds4/mission/apollo/v1/PDS4_APOLLO_1G00_1000.xsd">
'''
    pattern = '<Product_Observational.*'
    xml_data = re.sub(pattern, repl, xml_data, count=0, flags=0)

    base = os.path.splitext(csv_filename)[0]
    xml_filename = base + '.xml'
    print('File to write ', xml_filename)

    # write (with correct line endings)
    xml_file = open(xml_filename, 'w', newline='\r\n')
    xml_file.write(xml_data)


def make_dir_lower(top_level_dir,network,data_stream_type,station,year,julday):
    # makes the directory if one is not found
    # /Users/cnunn/mars_data/Archive_example_from_Renee/PDS_output/data/xb/continuous_waveform/elyse/2019/273/

    network = network.lower()
    data_stream_type = data_stream_type.lower()
    station = station.lower()
    # year and julday should be strings 

    directory = find_dir_lower(top_level_dir,network,data_stream_type,station,
      year,julday)
    if os.path.exists(directory):
        return directory
    # check that the overall directory exists
    elif not os.path.exists(top_level_dir):
        msg = ("The directory {} doesn't exist".format(top_level_dir))
        raise IOError(msg)
    else:
        directory = os.path.join(top_level_dir, network)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, data_stream_type)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, station)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, year)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, julday)
        if not os.path.exists(directory):
            os.makedirs(directory)

        return directory


def make_metatdata_dir_lower(top_level_dir,network,data_stream_type,sta):
    # makes the directory if one is not found
    # /Users/cnunn/mars_data/Archive_example_from_Renee/PDS_output/data/xb/continuous_waveform/elyse/2019/273/

    network = network.lower()
    data_stream_type = data_stream_type.lower()
    sta = sta.lower()

    directory = find_metadata_dir_lower(top_level_dir,network,data_stream_type,sta)
    if os.path.exists(directory):
        return directory
    # check that the overall directory exists
    elif not os.path.exists(top_level_dir):
        msg = ("The directory {} doesn't exist".format(top_level_dir))
        raise IOError(msg)
    else:
        directory = os.path.join(top_level_dir, network)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, data_stream_type)
        if not os.path.exists(directory):
            os.makedirs(directory)
        directory = os.path.join(directory, sta)
        if not os.path.exists(directory):
            os.makedirs(directory)

        return directory

def find_dir(top_level_dir,network,data_stream_type,station,year,julday):
    # find and return the correct directory
    # /Users/cnunn/mars_data/Archive_example_from_Renee/PDS_output/data/xb/continuous_waveform/elyse/2019/273/
    return os.path.join(top_level_dir,network,data_stream_type,station,year,julday)

def find_dir_lower(top_level_dir,network,data_stream_type,station,year,julday):
    # find and return the correct directory in lower case
    network = network.lower()
    data_stream_type = data_stream_type.lower()
    station = station.lower()
    # /Users/cnunn/mars_data/Archive_example_from_Renee/PDS_output/data/xb/continuous_waveform/elyse/2019/273/
    return os.path.join(top_level_dir,network,data_stream_type,station,year,julday)

def find_metadata_dir_lower(top_level_dir,network,data_stream_type,sta):
    # find and return the correct directory in lower case
    network = network.lower()
    data_stream_type = data_stream_type.lower()
    sta = sta.lower()
    return os.path.join(top_level_dir,network,data_stream_type,sta)
# 
def archive(release_number):
    # make the metadata files - https://service.iris.edu/fdsnws/station/docs/1/builder/

    log_text = 'Start: ', datetime.now()
    logging.info(log_text)
    print(log_text)

    # make a stationxml file which contains ATT 
    # make_att_metadata(input_metadata,updated_metadata)
    # print('make metadata file (for geocsv processing), using latest file')
    # make_metadata(updated_metadata,metadata_file)

    # print('temporary - put these back')

    # 
    # print('make_dataless_labels')
    # make_dataless_labels(release_number)
    # 
    # print('make_stationxml_labels')
    # make_stationxml_labels(release_number)
    # #
    # # # print('copy_seed_files')
    # # # copy_seed_files(release_number)
    # # 
    print('make_seed_labels')
    make_seed_labels(release_number)
    # #
    # print('make_geocsv_files')
    # make_geocsv_files(release_number)
    print('make_geocsv_labels')
    make_geocsv_labels(release_number)

    log_text = 'End: ', datetime.now()
    logging.info(log_text)
    print(log_text)

def make_geocsv_labels(release_number):

    # start empty geocsv collection  
    geocsv_col = OrderedDict()



    for (geocsv_filepath, geocsv_basename) in find_files(OUTPUT_DIR, '*.a.csv'):
    # for (geocsv_filepath, geocsv_basename) in find_files(OUTPUT_DIR, 'xa.s12..att.1976.061.0.a.csv'):

        geocsv_filename = os.path.join(geocsv_filepath, geocsv_basename)

        print(geocsv_filename)

        path_list = geocsv_filepath.split(os.sep)

        julday = path_list[-1]
        year = path_list[-2]
        station = path_list[-3]
        data_stream_type = path_list[-4]
        network = path_list[-5]
        # 
        # temporary 
        # stop processing if data_stream_type != 'continuous_waveform'
        if data_stream_type.upper() != 'CONTINUOUS_WAVEFORM' :
            print(data_stream_type)
            break
            
        output_dir = make_dir_lower(OUTPUT_DIR,network=network,data_stream_type=data_stream_type,station=station,year=year,julday=julday)
        # # 
        # if geocsv_filename=='/Users/cnunn/lunar_data/PDART_PDS_output/xa/continuous_waveform/s12/1976/066/xa.s12.01.mhz.1976.066.1.1.a.csv':
        make_geocsv_xml(output_dir,geocsv_basename,release_number,geocsv_col)

            # print(geocsv_col)

        if TEST:
            # test a single file
            break

    write_collection_geocsv_file(geocsv_col)




def make_geocsv_files(release_number):
    '''
    Make the geocsv files using mseed2ascii. 

    The code is available from IRIS - https://github.com/iris-edu/mseed2ascii

    Note that the ATT files require at least 3 decimal places. 
    Use a version of mseed2ascii which is 2.7 or greater. 

    22c22
    < #define VERSION "2.6"
    ---
    > #define VERSION "2.7DEV"
    599c599
    <                                      "%-10.8g  ", *(float *)sptr);
    ---
    >                                      "%-10.9g  ", *(float *)sptr);
    602c602
    <                                      "%.8g", *(float *)sptr);
    ---
    >                                      "%.9g", *(float *)sptr);
    608c608
    <                                      "%-10.10g  ", *(double *)sptr);
    ---
    >                                      "%-10.17g  ", *(double *)sptr);
    611c611
    <                                      "%.10g", *(double *)sptr);
    ---
    >                                      "%.17g", *(double *)sptr);
    680c680
    <         outsize += snprintf (outbuffer + outsize, sizeof(outbuffer) - outsize, "%s%s%s %.8g\n",
    ---
    >         outsize += snprintf (outbuffer + outsize, sizeof(outbuffer) - outsize, "%s%s%s %.9g\n",
    684c684
    <         outsize += snprintf (outbuffer + outsize, sizeof(outbuffer) - outsize, "%s%s%s %.10g\n",
    ---
    >         outsize += snprintf (outbuffer + outsize, sizeof(outbuffer) - outsize, "%s%s%s %.17g\n",
    688a689,690
    >       printf("%d", 1000);

    Make the file in the source directory and then make a symbolic link to it here: 
    ln -s /usr/local/mseed2ascii-chad2021/mseed2ascii mseed2ascii 

    > 

    '''

    for (seed_filepath, seed_basename) in find_files(OUTPUT_DIR, '*.mseed'):
    # for (seed_filepath, seed_basename) in find_files(OUTPUT_DIR, '*.mseed'):
        print(seed_filepath, seed_basename)
        # seed_basename = 'xb.elyse.02.bhw.2019.307.7.mseed'
        # seed_filepath='/Users/cnunn/mars_data/Archive/output_for_PDS/data_delivery_release_4/xb/continuous_waveform/elyse/2019/307'
        seed_filename = os.path.join(seed_filepath, seed_basename)

        base = os.path.splitext(seed_filename)[0]
        geocsv_tmp = '{}.a.tmp'.format(base)
        geocsv_filename = '{}.a.csv'.format(base)

        path_list = seed_filepath.split(os.sep)
        
        # /Users/cnunn/lunar_data/PDART_PDS_output/xa/continuous_waveform/s17/1976/067
        # xa.s17.9.mhz.1976.067.1.mseed
        julday = path_list[-1]
        year = path_list[-2]
        station = path_list[-3]
        data_stream_type = path_list[-4]
        network = path_list[-5]

        # print('temporary')
        # if 'afr' in seed_filename:
        #     print('afr found')
        #     continue
        # if 'att' in seed_filename:
        #     print('att found')
        #     continue



        # stop processing if data_stream_type != 'continuous_waveform'
        if data_stream_type.upper() != 'CONTINUOUS_WAVEFORM' :
            print(data_stream_type)
            break
        
        # to put extra information into the GeoCSV file, use something like this
        # p = Popen(['mseed2ascii','-G','-f','2', '-m' , metadata_file, seed_filename, '-o',geocsv_tmp, '-dr'], stdout=PIPE, stderr=PIPE)

        # just use the basic information in the header of the GeoCSV file - the extra information
        # is in the metadata file (note that we need to use different metadata files for the )
        p = Popen(['mseed2ascii','-G','-f','2', seed_filename, '-o',geocsv_tmp, '-dr'], stdout=PIPE, stderr=PIPE)
        # print('mseed2ascii','-G','-f','2', seed_filename, '-o',geocsv_tmp, '-dr', 'stdout=PIPE', 'stderr=PIPE')

        stdout, stderr = p.communicate()
        stderr_str = stderr.decode('utf-8') 
        
        # some examples for testing
        # stderr_str = 'Reported sample rate different than derived rate (20 versus 20)'
        # stderr_str = 'Reported sample rate different than derived rate (30 versus 20)'
        # stderr_str = 'Reported sample rate different than derived rate (1 versus 0.98971)'
        # stderr_str = 'Wrote 43200 samples for' 

        # Do some error handling
        # if we find something like 'Wrote 43200 samples for', then it's good 
        pattern = re.compile(r'^Wrote*')
        wrote = pattern.search(str(stderr_str))
        if wrote is None:
            # We can also check when it says the reported rate is different from 
            # the derived rate. If it's only a tiny change, we will use the 
            # reported rate
            pattern = re.compile(r'\(.*?versus.*?\)')
            versus = pattern.search(str(stderr_str))
            if versus:
                print('Ignoring error message and continuing: {}\n{}'.format(seed_filename, stderr_str))
                    # first_found = versus.group(0)
                    # array = re.findall(r'\d*[.,]?\d*', first_found) 
                    # number1 = None
                    # number2 = None
                    # for n in array:
                    #     try:
                    #         number = float(n)
                    #         if number1 is None:
                    #             number1 = number
                    #         else:
                    #             number2 = number
                    #     except:
                    #         pass
                    # if abs(100*((number1 - number2)/number1)) > 0.5:
                    #     print(seed_filename)
                    #     print(geocsv_filename)
                    #     print('NOT !!!!Exiting: {}'.format(stderr_str))
                    # else:
                    #     print(seed_filename)
                    #     print('Difference is less than 0.5%, so ignoring: {}'.format(stderr_str))
            else: 
                # the error is unknown, so print it and exit
                print(seed_filename)
                print(stderr_str)
                exit()



        # rewrite the line endings 
        tmp_file = open(geocsv_tmp, 'r')
        csv_file = open(geocsv_filename, 'w', newline='\r\n')
        
        for line in tmp_file:
            # if 'field_unit' in line:
                # also write out the geodatic_datum line
                # csv_file.write('# geodetic_datum: Mars 2000 planetocentric, MOLA geoid\n')
            csv_file.write(line)
        tmp_file.close()
        csv_file.close()
        os.remove(geocsv_tmp)
        print('Wrote file ', geocsv_filename)
        
        if TEST:
            # test a single file
            break

def make_dataless_labels(release_number):

    # make a collection of dataless files
    dataless_col = OrderedDict()

    for (meta_filepath, meta_basename) in find_metadata_output_files(OUTPUT_DIR, 'dataless*', case='lower'):
        print(meta_filepath, meta_basename)
        
        meta_filename = os.path.join(meta_filepath, meta_basename)
        print(meta_filename)

        path_list = meta_filepath.split(os.sep)
        # Found filepath: /Users/cnunn/mars_data/Archive/input_from_IPGP/data_delivery_release_4/XB/CONTINUOUS_WAVEFORM/ELYHK/2019/307
        # ['', 'Users', 'cnunn', 'mars_data', 'Archive', 'input_from_IPGP', 'data_delivery_release_4', 'XB', 'CONTINUOUS_WAVEFORM', 'ELYHK', '2019', '307']

        data_stream_type = path_list[-1]
        network = path_list[-2]

        print(data_stream_type, network)
        # 
        # temporary 
        # stop processing if data_stream_type != 'continuous_waveform'
        if data_stream_type.upper() != 'METADATA' :
            print(data_stream_type)
            break


        make_dataless_xml(meta_filepath,meta_basename,release_number,dataless_col)

        if TEST:
            # test a single file
            break

    write_collection_dataless_file(dataless_col)

def make_stationxml_labels(release_number):

    # new collection
    stationxml_col = OrderedDict()

    # /Users/cnunn/lunar_data/PDART_PDS_output/xa/metadata/stationxml.xa.2021.060.sxml
    for (meta_filepath, meta_basename) in find_metadata_output_files(OUTPUT_DIR, 'stationxml*.sxml', case='lower'):
        print(meta_filepath, meta_basename)
        meta_filename = os.path.join(meta_filepath, meta_basename)
        # print(seed_filename)

        path_list = meta_filepath.split(os.sep)
        # Found filepath: /Users/cnunn/mars_data/Archive/input_from_IPGP/data_delivery_release_4/XB/CONTINUOUS_WAVEFORM/ELYHK/2019/307
        # ['', 'Users', 'cnunn', 'mars_data', 'Archive', 'input_from_IPGP', 'data_delivery_release_4', 'XB', 'CONTINUOUS_WAVEFORM', 'ELYHK', '2019', '307']

        data_stream_type = path_list[-1]
        network = path_list[-2]
        # 
        # temporary 
        # stop processing if data_stream_type != 'continuous_waveform'
        if data_stream_type.upper() != 'METADATA' :
            print(data_stream_type)
            break

        make_stationxml_xml(meta_filepath,meta_basename,release_number,stationxml_col)

        if TEST:
            # test a single file
            break

    write_collection_stationxml_file(stationxml_col)



def make_seed_labels(release_number):

    # start empty seed collection
    seed_col = OrderedDict()

    for (seed_filepath, seed_basename) in find_files(OUTPUT_DIR, '*.mseed'):
        print(seed_filepath, seed_basename)

        seed_filename = os.path.join(seed_filepath, seed_basename)
        print(seed_filename)

        path_list = seed_filepath.split(os.sep)
        # /Users/cnunn/lunar_data/PDART_PDS_output/xa/continuous_waveform/s17/1976/067 

        # julday = path_list[-1]
        # year = path_list[-2]
        # station = path_list[-3]
        data_stream_type = path_list[-4]
        # network = path_list[-5]

        # 
        # temporary 
        # stop processing if data_stream_type != 'continuous_waveform'
        if data_stream_type.upper() != 'CONTINUOUS_WAVEFORM' :
            print(data_stream_type)
            break

        make_seed_xml(seed_filepath,seed_basename,release_number,seed_col)

        if TEST:
            # test a single file
            break


    write_collection_seed_file(seed_col)


        

        # input_filename = os.path.join(input_filepath, input_basename)
        # 
        # path_list = input_filepath.split(os.sep)
        # # Found filepath: /Users/cnunn/mars_data/Archive/input_from_IPGP/data_delivery_release_4/XB/CONTINUOUS_WAVEFORM/ELYHK/2019/307
        # # ['', 'Users', 'cnunn', 'mars_data', 'Archive', 'input_from_IPGP', 'data_delivery_release_4', 'XB', 'CONTINUOUS_WAVEFORM', 'ELYHK', '2019', '307']
        # 
        # data_stream_type = path_list[-1]
        # network = path_list[-2]
        # # 
        # # temporary 
        # # stop processing if data_stream_type != 'continuous_waveform'
        # if data_stream_type.upper() != 'METADATA' :
        #     print(data_stream_type)
        #     break
        # 
        # output_dir = make_metatdata_dir_lower(OUTPUT_DIR,network=network,data_stream_type=data_stream_type)
        # 
        # # dataless.xb.elyhk.2020.156.1.seed
        # filetype, net, sta, year, jul, _, ext = input_basename.split('.')
        # meta_basename = '{}.{}.{}.{}.{}.{}'.format(filetype, net, sta, year, jul, ext)
        # meta_basename = meta_basename.lower()
        # meta_filename = os.path.join(output_dir, meta_basename)
        # print(input_filename,meta_basename)
        # 
        # # base = os.path.splitext(seed_filename)[0]
        # shutil.copy2(input_filename,meta_filename)



def copy_seed_files(release_number):
    

    # for (input_filepath, input_basename) in find_initial_files(INPUT_DIR, '*S11*1969.202.*.MSEED', case='upper'):
    # for (input_filepath, input_basename) in find_initial_files(INPUT_DIR, '*1976.070*.MSEED', case='upper'):
    for (input_filepath, input_basename) in find_initial_files(INPUT_DIR, '*.MSEED', case='upper'):

        input_filename = os.path.join(input_filepath, input_basename)
        print(input_filename)

        path_list = input_filepath.split(os.sep)

        
        # Found filepath: /Users/cnunn/mars_data/Archive/input_from_IPGP/data_delivery_release_4/XB/CONTINUOUS_WAVEFORM/ELYHK/2019/307
        # ['', 'Users', 'cnunn', 'mars_data', 'Archive', 'input_from_IPGP', 'data_delivery_release_4', 'XB', 'CONTINUOUS_WAVEFORM', 'ELYHK', '2019', '307']
        # print(path_list)

        julday = path_list[-1]
        year = path_list[-2]
        station = path_list[-3]
        data_stream_type = 'continuous_waveform'
        network = 'XA'



        # # temporary! just to test a few
        # if year == 1976 and julday in ['061', '062', '063']:
        #     pass
        # else:
        #     continue

        # 
        # temporary 
        # stop processing if data_stream_type != 'continuous_waveform'
        # if data_stream_type.upper() != 'CONTINUOUS_WAVEFORM' :
        #     print(data_stream_type)
        #     break

        output_dir = make_dir_lower(OUTPUT_DIR,network=network,data_stream_type=data_stream_type,station=station,year=year,julday=julday)
        # print(output_dir)
        seed_basename = input_basename.lower()
        seed_filename = os.path.join(output_dir, seed_basename)


        # base = os.path.splitext(seed_filename)[0]
        shutil.copy2(input_filename,seed_filename)
        print(seed_filename)

        if TEST:
            # test a single file
            break

# does the file exist in the dictionary?

def update_collection_seed_file(seed_col,network,station,location,channel,year,julday,rev,version_id):
    # P,urn:nasa:pds:apollo_pse:data_seed:xa.s12.01.afr.1976.066.1.1::1.0
    root = '{}:{}.{}.{}.{}.{}.{}.{}::'.format(LOGICAL_IDENTIFIER_ROOT,network,station,location,channel,year,julday,rev)
    seed_col.update({root: version_id})
    return seed_col

def update_collection_geocsv_file(geocsv_col,network,station,location,channel,year,julday,rev,version_id):
    # P,urn:nasa:pds:apollo_pse:data_seed:xb.elyse.85.llz.2019.273.6::1.0
    # P,urn:nasa:pds:apollo_pse:data_table:xa.s12.01.mhz.1976.066.1.1.a::1.0
    root = '{}:{}.{}.{}.{}.{}.{}.{}.a::'.format(LOGICAL_IDENTIFIER_ROOT_GEOCSV,network,station,location,channel,year,julday,rev)
    geocsv_col.update({root: version_id})
    return geocsv_col

def update_collection_dataless_file(dataless_col,network,rev,version_id):
    # collection_data_seed_inventory.csv:
# P,urn:nasa:pds:apollo_pse:data_seed:dataless.xa.2021.060::1.0
# P,urn:nasa:pds:insight_seis:data_seed:dataless.xb.elys0.2019.144::1.1
    root = '{}:dataless.{}.{}::'.format(LOGICAL_IDENTIFIER_ROOT,network,rev)

    dataless_col.update({root: version_id})
    return dataless_col

def update_collection_stationxml_file(stationxml_col,network,rev,version_id):
    # P,urn:nasa:pds:apollo_pse:data_table:stationxml.xa.2021.060::1.0
    root = '{}:stationxml.{}.{}::'.format(LOGICAL_IDENTIFIER_ROOT_GEOCSV,network,rev)
    stationxml_col.update({root: version_id})
    return stationxml_col


def write_collection_seed_file(seed_col):
    with open(os.path.join(PREP_DIR, 'collection_data_seed_inventory.csv.tmp'),'w', newline='\r\n') as fp:
        for key in seed_col:
            # P,urn:nasa:pds:apollo_pse:data_seed:xb.elyse.85.llz.2019.273.6::1.0
            line1 = 'P,{}{}\n'.format(key, seed_col[key])
            fp.write('P,{}{}\n'.format(key, seed_col[key]))
    return seed_col

def write_collection_geocsv_file(geocsv_col):
    with open(os.path.join(PREP_DIR, 'collection_data_table_inventory.csv.tmp'),'w', newline='\r\n') as fp:
        for key in geocsv_col:
            # P,urn:nasa:pds:apollo_pse:data_table:xb.elyse.85.llz.2019.273.6.a::1.0
            line1 = 'P,{}{}\n'.format(key, geocsv_col[key])
            fp.write('P,{}{}\n'.format(key, geocsv_col[key]))
    return geocsv_col

def write_collection_dataless_file(dataless_col):
    with open(os.path.join(PREP_DIR, 'collection_data_seed_inventory_dataless.csv.tmp'),'w', newline='\r\n') as fp:
        for key in dataless_col:
            # P,urn:nasa:pds:apollo_pse:data_seed:dataless.xb.elyse.2020.062::1.0
            line1 = 'P,{}{}\n'.format(key, dataless_col[key])
            fp.write('P,{}{}\n'.format(key, dataless_col[key]))
    return dataless_col

def write_collection_stationxml_file(stationxml_col):
    with open(os.path.join(PREP_DIR, 'collection_data_table_inventory_stationxml.csv.tmp'),'w', newline='\r\n') as fp:
        for key in stationxml_col:
            # P,urn:nasa:pds:apollo_pse:data_seed:dataless.xb.elyse.2020.062::1.0
            line1 = 'P,{}{}\n'.format(key, stationxml_col[key])
            fp.write('P,{}{}\n'.format(key, stationxml_col[key]))
    return stationxml_col




# TODO Find an easier way to sort the results, so that you go through systematically - 
# At the moment need one find_files for each type of directory structure 
def find_files(directory, pattern, unmatch_pattern=None, network='xa', waveform_type='continuous_waveform', case='lower'):
    # find the files which match a pattern
    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()
    # print(os.path.join(directory,network,waveform_type))
    for station in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
        for year in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
            for julday in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,year))):
                for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,year,julday))):
                    # pass
                    if fnmatch.fnmatch(basename, pattern):
                        if unmatch_pattern is not None:
                            if not fnmatch.fnmatch(basename, unmatch_pattern):
                                # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
                                yield os.path.join(directory,network,waveform_type,station,year,julday), basename
                        else:
                            yield os.path.join(directory,network,waveform_type,station,year,julday), basename

def find_initial_files(directory, pattern, unmatch_pattern=None, network='xa', waveform_type='continuous_waveform', case='lower'):
    # find the files which match a pattern
    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()
    # print(os.path.join(directory,network,waveform_type))
    for station in sorted(os.listdir(directory)):
        for year in sorted(os.listdir(os.path.join(directory,station))):
            for julday in sorted(os.listdir(os.path.join(directory,station,year))):
                for basename in sorted(os.listdir(os.path.join(directory,station,year,julday))):
                    # pass
                    if fnmatch.fnmatch(basename, pattern):
                        if unmatch_pattern is not None:
                            if not fnmatch.fnmatch(basename, unmatch_pattern):
                                # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
                                yield os.path.join(directory,station,year,julday), basename
                        else:
                            yield os.path.join(directory,station,year,julday), basename

def find_metadata_files(directory, pattern, network='xa', waveform_type='metadata', case='upper'):
    # find the files which match a pattern

    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()
    # for station in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
    #     for year in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
    #         for julday in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,year))):
    # for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
    #                 # pass
    #     if fnmatch.fnmatch(basename, pattern):
    #         # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
    #         yield os.path.join(directory,network,waveform_type), basename

    # print(directory,network,waveform_typ)
    # print(os.path.join(directory,network,waveform_type))
    # print(os.path.join(directory,network,waveform_type))
    # print(sorted(os.listdir(os.path.join(directory,network,waveform_type))))
    for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type))):

        if fnmatch.fnmatch(basename, pattern):
            # print(os.path.join(directory,network,waveform_type,station), basename)
            yield os.path.join(directory,network,waveform_type), basename

def find_metadata_output_files(directory, pattern, network='xa', waveform_type='metadata', case='lower'):
    # find the files which match a pattern

    print('directory', directory )
    print('pattern ', pattern)

    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()

    for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
        if fnmatch.fnmatch(basename, pattern):
            # print(os.path.join(directory,network,waveform_type,station), basename)
            yield os.path.join(directory,network,waveform_type), basename

def find_event_files(directory, pattern, unmatch_pattern=None, network='xa', waveform_type='event_waveform', case='lower'):
    # find the files which match a pattern
# /Users/cnunn/mars_data/Archive/input_from_PDS/data_delivery_release_5a/xb/event_waveform/elyse/mqs2019fddj

    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()
    # print(os.path.join(directory,network,waveform_type))
    # for station in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
    #     for year in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
    #         for julday in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,year))):
    for station in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
        for event_id in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
            for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,event_id))):
                    # pass
                if fnmatch.fnmatch(basename, pattern):
                    if unmatch_pattern is not None:
                        if not fnmatch.fnmatch(basename, unmatch_pattern):
                            # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
                            yield os.path.join(directory,network,waveform_type,station,event_id), basename
                    else:
                        # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
                        yield os.path.join(directory,network,waveform_type,station,event_id), basename
                    
                    


def find_metadata2_files(directory, pattern, network='xa', waveform_type='metadata', case='lower'):
    # find the files which match a pattern

    if case == 'upper':
        waveform_type = waveform_type.upper()
        network = network.upper()
    # print(os.path.join(directory,network,waveform_type))
    for station in sorted(os.listdir(os.path.join(directory,network,waveform_type))):
    #     for year in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
    #         for julday in sorted(os.listdir(os.path.join(directory,network,waveform_type,station,year))):
        for basename in sorted(os.listdir(os.path.join(directory,network,waveform_type,station))):
                        # pass
            if fnmatch.fnmatch(basename, pattern):
                # print(os.path.join(directory,network,waveform_type,station,year,julday, basename))
                yield os.path.join(directory,network,waveform_type,station), basename

def make_att_metadata(input_metadata,updated_metadata):
    # Add AFR and ATT to the metadata files 
    inv = read_inventory(input_metadata)
    #
    network = inv[0]
    for station in network:
        # print(station)
        channel_att = None
        for channel in station:
            if channel.code == 'MHZ':
                # from obspy.core.inventory.channel import Channel
                # channel1 = Channel()
                # print(channel1.depth)

                channel_att = channel.copy()
                channel_att.code = 'ATT'
                channel_att.sensor.description = 'Time the signal was recorded on Earth at the ground station. Time in seconds since the beginning of the epoch (1 January 1970). The timing can be negative for early in the mission.'
                # no response
                channel_att.response = None
                # channel startdate/enddate equals station startdate/enddate
                channel_att.start_date = station.start_date
                channel_att.end_date = station.end_date
                channel_att.comments = []
                channel_att.location_code = ''
                channel_att.azimuth=None
                channel_att.dip=None
                channel_att.sample_rate=1/(DELTA*4)

                # print(channel_att)
                # exit()
                # only need one per station
                break 


        station.channels.append(channel_att)

        print(inv)

    inv.write(updated_metadata,format='STATIONXML')
    print('Metadata file written ', updated_metadata)


# obspy.io.xseed.scripts.xseed2dataless
# https://docs.obspy.org/packages/autogen/obspy.io.xseed.scripts.xseed2dataless.html


def make_metadata(input_metadata,metadata_file):
    # this is unnecessary - no longer putting the extra meta information
    # into the GeoCSV files 


    # Can also check the metadata required by the geocsv files (when it's at IRIS)
    # https://service.iris.edu/fdsnws/station/docs/1/builder/

    # input metadata can be either dataless or stationxml
    # very large files take several seconds to process
    # input_metadata = '/Users/cnunn/mars_data/Archive/output_for_PDS/data_delivery_release_7/xb/metadata/elyse/dataless.xb.elyse.2020.345.seed'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.xb.elyse.2020.345.csv'
    # 
    # input_metadata = '/Users/cnunn/mars_data/Archive/output_for_PDS/data_delivery_release_7/xb/metadata/elyse/dataless.xb.elyse.2020.345.seed'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.xb.elyse.2020.345.csv'
    # 
    # # original lunar dataless file 
    # input_metadata = '/Users/cnunn/lunar_data/PDART_PDS_output/xa/metadata/XA.1969-1977_original.dataless'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.XA.1969-1977_original.csv'
    # 
    # # improved lunar stationxml (may still be some errors)
    # input_metadata = '/Users/cnunn/lunar_data/Dataless_SEED/XA.1969-1977_updated_2019.xml'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.XA.1969-1977_updated_2019.csv'

    

    with open(metadata_file, 'w', newline='\r\n') as fm:
        inv = read_inventory(input_metadata)
        print(inv)
        header = '#Network | Station | Location | Channel | Latitude | Longitude | Elevation | Depth | Azimuth | Dip | SensorDescription | Scale | ScaleFreq | ScaleUnits | SampleRate | StartTime | EndTime'
        fm.write(header)
        network = inv[0]
        for station in network:
            # print(station)
            for channel in station:
                # print(channel)
                # print(channel.response)
                description = channel.sensor.description
                if description is None:
                    description = channel.sensor.type

                if channel.response is not None:
                    scale = '{:.5E}'.format(channel.response.instrument_sensitivity.value)
                    scale_freq = channel.response.instrument_sensitivity.frequency
                    scale_units = channel.response.instrument_sensitivity.input_units.lower()
                else:
                    scale = ''
                    scale_freq = ''
                    scale_units = ''

                line = '{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}\n'.format(
                    network.code,station.code,channel.location_code,
                    channel.code,station.latitude,station.longitude,
                    station.elevation,channel.depth,channel.azimuth,
                    channel.dip,description,
                    scale, scale_freq, scale_units,
                    channel.sample_rate,
                    channel.start_date.strftime('%Y-%m-%dT%H:%M:%S'), 
                    channel.end_date.strftime('%Y-%m-%dT%H:%M:%S'))

                fm.write(line)
        print('Metadata file written ', metadata_file)

    

# InSight example
# Should look like this:
# XB|ELYSE|00|HHU|4.502384|135.623447|-2613.4|-0.1|135.1|-29.4|VBB Velocity SCI mode|8.62622E10|0.1|m/s|100.0|2019-02-10T14:16:51|2019-02-11T21:46:17
# Does look like this (the same except for the shorthand scientific notation in the previous one)
# XB|ELYSE|00|HHU|4.502384|135.623447|-2613.4|-0.1|135.1|-29.4|VBB Velocity SCI mode|86262200000.0|0.1|m/s|100.0|2019-02-10T14:16:51|2019-02-11T21:46:17

# Apollo examples 
# original dataless file 
# # XA|S12||MH1|-3.04|-23.42|123456.0|0.0|180.0|0.0|Apollo PSE Alsep Seismometer|3579160000.0|1.0|m|6.625|1969-11-19T14:23:00|1974-10-16T14:02:00
# The corrected one (as a stationxml)
# # XA|S12||MH1|-3.01084|-23.42456|-1424.0|0.0|180.0|0.0|Apollo PSE Alsep Seismometer|3579160000.0|1.0|m|6.625|1969-11-19T14:23:00|1974-10-16T14:02:00

        

if __name__ == "__main__":

    



    # TODO check all files have correct line endings 

    # copy the metadata file
    # cp /Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0.xml /Users/cnunn/lunar_data/PDART_PDS_OUTPUT/data/xa/metadata/stationxml.xa.0.sxml

    # make the dataless file 
    # cd /Users/cnunn/lunar_data/PDART_PDS_OUTPUT/data/xa/metadata
    # java -jar /Users/cnunn/Applications/stationxml-seed-converter-2.1.0.jar --input stationxml.xa.0.sxml --output dataless.xa.0.seed

    # release '1' 
    release_number='1'
    INPUT_DIR = '/Users/cnunn/lunar_data/PDART/'
    OUTPUT_DIR = '/Users/cnunn/lunar_data/pdart_pds_output/data/'
    PREP_DIR = '/Users/cnunn/lunar_data/PDART_PDS_PREP'

    # input_metadata = '/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0X.xml'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.XA.1969-1977_updated_2019.csv'
    # updated_metadata = '/Users/cnunn/lunar_data/PDART_METADATA/XA.1969-1977.0Y.xml'
    # metadata_file = '/Users/cnunn/python_packages/apollo_archive/files/metadata.XA.1969-1977.0.csv'

    archive(release_number)





    # stream = read('/Users/cnunn/lunar_data/PDART/1970/XA/S12/ATT/XA.S12.02.ATT.1970.050.gz')
    # stream = read('/Users/cnunn/lunar_data/PDART/1976/XA/S12/ATT/XA.S12.15.ATT.1976.055.gz')
    # 
    # stream = read('/Users/cnunn/lunar_data/PDART_PDS_output/xa/continuous_waveform/s17/1976/066/xa.s17.06.mhz.1976.066.1.1.mseed')
    # 
    # print(stream[0].data[42])
    # # print(UTCDateTime(stream[0].data[0]))
    # print(stream[0])
    # exit()


