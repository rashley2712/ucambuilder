#!/usr/bin/env python

"""
Some helper routines for processing form data sent by email
"""

import re
import exceptions
import trm.subs as subs

def ttime(tstr):
    """
    ttime translates the time string signifying arrival time into a float
    so that the messages can be ordered by arrival time.
    """

    days = (0,31,59,90,120,151,181,212,243,273,304,334)
    try:
        wday,mday,month,year,utc = tstr.split()[:5]
        (hour,min,sec) = utc.split(':')
        iy = int(year)
        if iy % 4 == 0:
            lday = 1
        else:
            lday = 0
        im = subs.month2int(month) - 1

        day = float(365*(iy - 2000)) + float(days[im]) + float(mday) - 1 + (int(hour) + int(min)/60. + int(sec)/3600.)/24.
        # correct for leap day of the current year
        if im > 1:
            day += 1.
        # correct for leap days of previous years
        for y in range(2000,iy):
            if y % 4 == 1:
                day += 1
    except:
        day = 0
            
    return day
        

def glue(str):
    """
    glue prevent HTML line breaking by sticking in non-breaking spaces
    """
    return '&nbsp;'.join(str.split())


def poskey(fields):
    """
    Returns a key out of the RA, Dec and target name for ordering targets.

    fields -- a dictionary with keys 'RA (J2000)' and 'Declination (J2000)' which return
              strings giving the corresponding values and a key 'Target name'.
    """
    ram = re.compile(r'^J?\s*(\d{1,2})[^0-9]+?(\d{1,2})[^0-9]+?(\d{1,2}(?:\.\d*)?)[^0-9]*$').match(fields['RA (J2000)'])
    rapart  = '%02d%02d%05.2f' % (int(ram.group(1)),int(ram.group(2)),float(ram.group(3)))
    decm = re.compile(r'^\s*(-?|\+?|\s?)\s*(\d{1,2})[^0-9]+?(\d{1,2})[^0-9]+?(\d{1,2}(\.\d*)?)[^0-9]*$').match(fields['Declination (J2000)'])
    decpart = '%s%02d%02d%04.1f%s' % (decm.group(1),int(decm.group(2)),int(decm.group(3)),float(decm.group(4)),fields['Target name'])
    return rapart + decpart

def ultparse(formdata):
    """
    ultparse parses the text part of the e-mails from the phase II forms for
    ultracam and ultraspec returning a dictionary of the field values.
    """

    lines = formdata.get_payload(decode=True).splitlines()
    # Delete everything before Programme ID
    for i in range(len(lines)):
        if lines[i].startswith('Programme ID:'):
            del lines[:i]
            break
    else:
        raise Exception('ultparse error: could not find Programme ID')

    # delete everthing after 'Click ..'
    for i in range(len(lines)):
        if lines[i].startswith('Click'):
            del lines[i-1:]
            break

    # combine comment strings
    first = True
    for i in range(len(lines)):
        if lines[i].startswith('Comments:') or lines[i].startswith('Comments specific to target:'):
            if i < len(lines)-1:
                for j in range(i+1,len(lines)):
                    if not lines[j].isspace() and lines[j] != '':
                        if first:
                            lines[i] += ' ' + lines[j]
                            first = False
                        else:
                            lines[i] += '<br> ' + lines[j]
                del lines[i+1:]
            break
    else:
        raise FormsError('ultparse error: could not find any comments')
    
    # list of keywords. Not all need to exist but these must cover all
    # possibilities.
    keywords = ( \
        'Programme ID','Your name','Your e-mail address','Telephone number',\
            'Alternative telephone number', 'Number of nights awarded', 'Yes', \
            'Your proposal (pdf or ps)', 'Comments', 'Target name', 'Priority', \
            'RA (J2000)', 'Declination (J2000)', 'Magnitude', 'Filters (e.g. ugr)', \
            'Desired exposure time (sec)','Maximum exposure time (sec)','OK to use u-band coadds',\
            'Total time on target (mins)','Minimum unbroken time on target (mins)','Must be photometric',\
            'Seeing (worst case)','Lunar phase (worst case)','Time critical (e.g. eclipses and satellites)',\
            'Maximum u-band exposure time (sec)','Finding chart(s) (see below)','Comments specific to target', \
            'Exposure time (sec) (S)', 'Total time on target (mins) (S)', 'Minimum unbroken time on target (mins) (S)', \
            'Standards', 'Slit angle', 'Slit', 'Grism', 'Filters (S)', 'Expected S-to-N', \
            'Minimum S-to-N', 'Filters (I)', 'Exposure time (sec) (I)', 'Total time on target (mins) (I)', \
            'Minimum unbroken time on target (mins) (I)')

    # stick into a dictionary, cope with a couple of special cases.
    field = {}
    j = 0
    for i in range(len(lines)):
        for key in keywords:
            if lines[j].startswith(key + ':'):
                if key.startswith('Comments'):
                    value = lines[j][lines[j].find(':')+1:]
                    field[key] = value
                else:
                    if j < len(lines)-1:
                        if key == 'Telephone number' or key == 'Alternative telephone number' or key == 'Your name':
                            field[key] = glue(lines[j+1].strip())
                        elif key == 'Total time on target (mins)' or key == 'Minimum unbroken time on target (mins)':
                            num = re.compile(r'(\d+(?:\.\d*)?)')
                            m = num.match(lines[j+1].strip())
                            if m:
                                field[key] = m.group(1)
                            else:
                                field[key] = '0'
                        elif key == 'Filters (e.g. ugr)':
                            field[key] = lines[j+1].strip().replace(',',', ')
                        else:
                            field[key] = lines[j+1].strip()
                        j += 1
                    else:
                        field[key] = ''
                break
        j += 1
        if j >= len(lines):
            break

    # return
    return field


def ultchk(fields):
    """
    checks various fields of the ultra-cam/-spec forms
    """
    progID = fields['Programme ID']
    rachk  = re.compile(r'^J?\s*\d{1,2}[^0-9]+?\d{1,2}[^0-9]+?\d{1,2}(\.\d*)?[^0-9]*$')
    decchk = re.compile(r'^\s*(\-|\+|\s)?\s*\d{1,2}[^0-9]+?\d{1,2}[^0-9]+?\d{1,2}(\.\d*)?[^0-9]*$')
    for (key,value) in fields.iteritems():
        if key == 'RA (J2000)' and not rachk.match(fields[key]):
            print 'Programme = ' + progID + ', target = ' + fields['Target name'] + ', RA = ' + value + ' has invalid format'
            print 'Will skip this entry.'
            return False
        if key == 'Declination (J2000)' and not decchk.match(fields[key]):
            print 'Programme = ' + progID + ', target = ' + fields['Target name'] + ', Declination = ' + value + ' has invalid format'
            print 'Will skip this entry.'
            return False
        if key == 'Total time on target (mins) (S)' or key == 'Total time on target (mins) (I)':
            try:
                f = float(value)
            except:
                # special case
                if value == '2 nights':
                    fields[key] = str(2.0*10.5*60)
                else:
                    print 'Programme = ' + progID + ', target = ' + fields['Target name'] + ', could not interpret total time = ' + value
                    fields[key] = str(0.0)
                    print 'Will set = 0.'         
    return True


def genparse(formdata, fields):
    """
    parses the text part of the e-mails returning a dictionary of the field values. This one is generic
    which means that it does not do much.
    formdata -- the e-mail returned by the formsbuilder form
    fields   -- the fields to look for. The first one should be the first field of the form which
                will define the start of it. Do not include the ending ':'
    """

    lines = formdata.get_payload(decode=True).splitlines()

    # Top and tail the message:
    for i in range(len(lines)):
        if lines[i].startswith(fields[0] + ':'):
            del lines[:i]
            break
    else:
        raise FormError('genparse error: could not find ' + fields[0])

    for i in range(len(lines)):
        if lines[i].startswith('Click <http'):
            del lines[i-1:]
            break

    # Define a dictionary
    fdict = {}
    j = 0
    for i in range(len(lines)):
        for key in fields:
            if lines[j].startswith(key + ':'):
                if j < len(lines)-1:
                    fdict[key] = lines[j+1].strip()
                    j += 1
                elif j == len(lines)-1:
                    fdict[key] = ''
                break
        j += 1
        if j >= len(lines): break

    return fdict

# Exception class
class FormError(exceptions.Exception):
    def __init__(self, value):
        self.value = value
            
    def __str__(self):
        return repr(self.value)
    
