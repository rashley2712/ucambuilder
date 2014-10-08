#!/usr/bin/env python

descr = \
"""
provides a class to access google spreadsheet data, worksheet by
worksheet. At the moment it only works on privately accessible 
documents so you will need access to the document through a gmail 
account. 
"""

import gdata.docs.service
import gdata.spreadsheet.service

class Ssheet(object):
    """
    Connects to a google spreadsheet to allow the worksheets
    to be downloaded

    The class is iterable and callable, thus the folliwng code 
    is valid:

    from trm import gdoc

    ....

    enter your gmail account details (email address and password)
    .....

    ssheet = gdoc.Ssheet(email, password, 'My Spread Sheet')

    # get a particular worksheet
    table = ssheet('A worksheet')

    # get all of them
    for name,table in ssheet:
       print 'Found worksheet',name
       ... do something with each table


    """
    
    def __init__(self, email, password, ssheet):
        """
        email    -- your gmail address
        password -- password
        ssheet   -- the spreadsheet to access

        Raises an Exception if the spreadsheet 
        can't be found. Invalid login details
        will also rasie Exceptions.
        """

        # Connect to the document list service
        self._gd_client = gdata.docs.service.DocsService(source="getsheet")
        self._gd_client.ClientLogin(email, password)

        # connect to the google spreadsheet service
        self._gs_client = gdata.spreadsheet.service.SpreadsheetsService()
        self._gs_client.ClientLogin(email, password)

        # track tokens (not sure this is needed)
        self._docs_auth_token   = self._gd_client.GetClientLoginToken()
        self._sheets_auth_token = self._gs_client.GetClientLoginToken()

        self._gd_client.SetClientLoginToken(self._sheets_auth_token)

        # get document feed
        self._gs_feed = self._gs_client.GetSpreadsheetsFeed()
        for entry in self._gs_feed.entry:
            if entry.title.text == ssheet:
                self._skey = entry.id.text.split('/')[-1]
                print 'Found spreadsheet =',ssheet
                break
        else:
            raise Exception('Ssheet.__init__: failed to find spreadsheet = ' + ssheet)

        # get its worksheets
        self._wsheets = self._gs_client.GetWorksheetsFeed(self._skey)

    def __iter__(self):
        """
        Iterator to deliver next worksheet.
        Returns (name,table) where name is the
        name of the worksheet amd table are
        its values as a list of lists.
        """

        for entry in self._wsheets.entry:
            wskey = entry.id.text.split('/')[-1]

            # get column names using a cell feed
            query = gdata.spreadsheet.service.CellQuery()
            query.max_row = '1'
            query.min_row = '1'
            cfeed = self._gs_client.GetCellsFeed(self._skey, wskey, query=query)
            cnames = [centry.content.text for centry in cfeed.entry]
            keys = [cname.replace(' ','').lower() for cname in cnames]
            
            # get list feed to access the table data
            lfeed = self._gs_client.GetListFeed(self._skey, wskey)

            table = [cnames,]
            for lentry in lfeed.entry:
                row = [lentry.custom[key].text for key in keys]
                table.append([ent.strip() if ent else '' for ent in row])

            yield (entry.title.text, table)

    def __call__(self, wsheet):
        """
        Returns the worksheet called wsheet.
        Raises an Exception if it can't be found.
        """

        for entry in self._wsheets.entry:
            if entry.title.text == wsheet:
                wskey = entry.id.text.split('/')[-1]

                # get column names using a cell feed
                query = gdata.spreadsheet.service.CellQuery()
                query.max_row = '1'
                query.min_row = '1'
                cfeed = self._gs_client.GetCellsFeed(self._skey, wskey, query=query)
                cnames = [centry.content.text for centry in cfeed.entry]

                # next is a bit rubbish and it would better to replace. Can column headings
                # 'A,B', 'ab' and 'A()B' really all be equivalent?
                keys = [cname.replace(' ','').replace('(','').replace(')','').replace(',','').lower() \
                            for cname in cnames]

                # get list feed to access the table data
                lfeed = self._gs_client.GetListFeed(self._skey, wskey)

                table = [cnames,]
                for lentry in lfeed.entry:
                    row = [lentry.custom[key].text for key in keys]
                    table.append([ent if ent else '' for ent in row])
                return table

        raise Exception('Ssheet.__call__: failed to find worksheet = ' + wsheet)            

