# ***********************************************************************
# ******************  CANADIAN ASTRONOMY DATA CENTRE  *******************
# *************  CENTRE CANADIEN DE DONNÉES ASTRONOMIQUES  **************
#
#  (c) 2025.                            (c) 2025.
#  Government of Canada                 Gouvernement du Canada
#  National Research Council            Conseil national de recherches
#  Ottawa, Canada, K1A 0R6              Ottawa, Canada, K1A 0R6
#  All rights reserved                  Tous droits réservés
#
#  NRC disclaims any warranties,        Le CNRC dénie toute garantie
#  expressed, implied, or               énoncée, implicite ou légale,
#  statutory, of any kind with          de quelque nature que ce
#  respect to the software,             soit, concernant le logiciel,
#  including without limitation         y compris sans restriction
#  any warranty of merchantability      toute garantie de valeur
#  or fitness for a particular          marchande ou de pertinence
#  purpose. NRC shall not be            pour un usage particulier.
#  liable in any event for any          Le CNRC ne pourra en aucun cas
#  damages, whether direct or           être tenu responsable de tout
#  indirect, special or general,        dommage, direct ou indirect,
#  consequential or incidental,         particulier ou général,
#  arising from the use of the          accessoire ou fortuit, résultant
#  software.  Neither the name          de l'utilisation du logiciel. Ni
#  of the National Research             le nom du Conseil National de
#  Council of Canada nor the            Recherches du Canada ni les noms
#  names of its contributors may        de ses  participants ne peuvent
#  be used to endorse or promote        être utilisés pour approuver ou
#  products derived from this           promouvoir les produits dérivés
#  software without specific prior      de ce logiciel sans autorisation
#  written permission.                  préalable et particulière
#                                       par écrit.
#
#  This file is part of the             Ce fichier fait partie du projet
#  OpenCADC project.                    OpenCADC.
#
#  OpenCADC is free software:           OpenCADC est un logiciel libre ;
#  you can redistribute it and/or       vous pouvez le redistribuer ou le
#  modify it under the terms of         modifier suivant les termes de
#  the GNU Affero General Public        la “GNU Affero General Public
#  License as published by the          License” telle que publiée
#  Free Software Foundation,            par la Free Software Foundation
#  either version 3 of the              : soit la version 3 de cette
#  License, or (at your option)         licence, soit (à votre gré)
#  any later version.                   toute version ultérieure.
#
#  OpenCADC is distributed in the       OpenCADC est distribué
#  hope that it will be useful,         dans l’espoir qu’il vous
#  but WITHOUT ANY WARRANTY;            sera utile, mais SANS AUCUNE
#  without even the implied             GARANTIE : sans même la garantie
#  warranty of MERCHANTABILITY          implicite de COMMERCIALISABILITÉ
#  or FITNESS FOR A PARTICULAR          ni d’ADÉQUATION À UN OBJECTIF
#  PURPOSE.  See the GNU Affero         PARTICULIER. Consultez la Licence
#  General Public License for           Générale Publique GNU Affero
#  more details.                        pour plus de détails.
#
#  You should have received             Vous devriez avoir reçu une
#  a copy of the GNU Affero             copie de la Licence Générale
#  General Public License along         Publique GNU Affero avec
#  with OpenCADC.  If not, see          OpenCADC ; si ce n’est
#  <http://www.gnu.org/licenses/>.      pas le cas, consultez :
#                                       <http://www.gnu.org/licenses/>.
#
#  $Revision: 4 $
#
# ***********************************************************************
#

"""
This module implements the ObsBlueprint mapping, as well as the workflow entry point that executes the workflow.
"""

import logging
from os.path import basename
from caom2 import ProductType
from caom2pipe.astro_composable import get_datetime_mjd
from caom2pipe import caom_composable as cc
from caom2pipe.manage_composable import make_datetime, StorageName


__all__ = [
    'mapping_factory',
    'JCMTSLName',
]


class JCMTSLName(StorageName):
    """Naming rules:
    - support mixed-case file name storage, and mixed-case obs id values
    - support uncompressed files in storage
    """

    JCMTSL_NAME_PATTERN = '*'

    def __init__(self, source_names):
        super().__init__(file_name=basename(source_names[0]), source_names=source_names)

    def is_valid(self):
        return True

    def product_type(self):
        result = ProductType.SCIENCE
        if 'um.err.' in self._file_name:
            result = ProductType.NOISE
        elif 'um.cov.' in self._file_name:
            result = ProductType.WEIGHT
        elif 'um.obj.' in self._file_name:
            result = ProductType.AUXILIARY
        return result

    def set_obs_id(self, **kwargs):
        self._obs_id = self._file_name.split('.')[0]


class JCMTSLBase(cc.TelescopeMapping2):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        # mapfits from https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/JCMTSL/scubacat.pdf
        bp.set('DerivedObservation.members', {})
        bp.set('Observation.algorithm.name', 'mapfits')
        bp.set('Observation.instrument.name', 'SCUBA')
        bp.set('Observation.intent', 'science')
        bp.set('Observation.target.name', '_get_target_name()')
        bp.set('Observation.target.type', 'field')
        bp.set('Observation.telescope.name', 'JCMT')
        bp.set('Observation.telescope.geoLocationX', -5461075.78)
        bp.set('Observation.telescope.geoLocationY', -2491090.15)
        bp.set('Observation.telescope.geoLocationZ', 2149569.763)

        bp.set('Plane.calibrationLevel', 3)

        bp.set('Artifact.productType', self._storage_name.product_type())
        bp.set('Artifact.releaseType', 'data')

        self._logger.debug('Done accumulate_bp.')

    def update(self):
        """Called to fill multiple CAOM model elements and/or attributes (an n:n relationship between TDM attributes
        and CAOM attributes).
        """
        return super().update()

    def _get_target_name(self, ext):
        bits = self._storage_name.file_name.split('_')
        return f'{bits[2]}_{bits[3]}'

    def _update_artifact(self, artifact):
        for part in artifact.parts.values():
            for chunk in part.chunks:
                # constructed from hard-coded values, no cutouts to support
                chunk.energy_axis = None
                chunk.time_axis = None


class JCMTSLSpatialSpectral(JCMTSLBase):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)

        bp.add_attribute('Plane.dataRelease', 'DATE-OBS')
        bp.add_attribute('Plane.metaRelease', 'DATE-OBS')
        bp.set('Plane.dataProductType', 'image')
        bp.set('Plane.provenance.name', 'JCMTSL')
        bp.set('Plane.provenance.lastExecuted', '_get_provenance_last_executed()')
        bp.set('Plane.provenance.producer', '_get_provenance_producer()')
        bp.set('Plane.provenance.reference', 'https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/scubalegacy/')
        bp.set('Plane.provenance.version', '_get_provenance_version()')

        bp.configure_position_axes((1, 2))
        # J2000 equatorial coordinates
        bp.set('Chunk.position.axis.axis1.cunit', 'deg')
        bp.set('Chunk.position.axis.axis2.cunit', 'deg')
        # bp.set('Chunk.position.axis.function.cd11', 0.6 / self._headers[0].get('NAXIS1') )
        # 1.2 is from the paper as the size of each field
        # bp.set('Chunk.position.axis.function.cd11', 1.2 / self._headers[0].get('NAXIS1') )
        bp.add_attribute('Chunk.position.axis.function.cd11', 'CDELT1' )
        bp.set('Chunk.position.axis.function.cd12', 0.0)
        bp.set('Chunk.position.axis.function.cd21', 0.0)
        # bp.set('Chunk.position.axis.function.cd22', 0.6 / self._headers[0].get('NAXIS2'))
        # bp.set('Chunk.position.axis.function.cd22', 1.2 / self._headers[0].get('NAXIS2'))
        bp.add_attribute('Chunk.position.axis.function.cd22', 'CDELT2')

        bp.configure_energy_axis(3)
        # from an existing JCMT product
        # specsys: BARYCENT
        # ssysobs: TOPOCENT
        # ssyssrc: TOPOCENT
        # ctype: WAVE
        # cunit: um
        bp.set('Chunk.energy.axis.axis.ctype', 'WAVE')
        bp.set('Chunk.energy.axis.axis.cunit', 'um')
        bp.set('Chunk.energy.specsys', 'BARYCENT')
        bp.set('Chunk.energy.ssysobs', 'TOPOCENT')
        bp.set('Chunk.energy.ssyssrc', 'TOPOCENT')
        bp.set('Chunk.energy.axis.range.start.pix', 0.5)
        bp.set('Chunk.energy.axis.range.end.pix', 1.5)
        self._logger.debug('Done accumulate_bp.')

    def _get_date_obs(self, ext):
        return get_datetime_mjd(self._headers[ext].get('DATE-OBS'))

    def _get_provenance_last_executed(self, ext):
        result = None
        history = self._headers[0].get('HISTORY')
        for line in history:
            if 'Executed on:' in line:
                result = make_datetime(line.split('Executed on:')[1].strip())
                break
        return result

    def _get_provenance_producer(self, ext):
        result = None
        origin = self._headers[0].get('ORIGIN')
        if origin:
            result = origin.split(': version')[0]
        return result

    def _get_provenance_version(self, ext):
        result = None
        origin = self._headers[0].get('ORIGIN')
        if origin:
            result = origin.split(': version')[1]
        return result


class JCMTSL850um(JCMTSLSpatialSpectral):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)
        # 6" from https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/JCMTSL/scubacat.pdf
        bp.set('Chunk.position.resolution', 6.0)
        # 19" from https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/JCMTSL/scubacat.pdf
        bp.set('Chunk.energy.axis.range.start.val', 840.5)
        bp.set('Chunk.energy.axis.range.end.val', 859.5)
        bp.set('Chunk.energy.bandpassName', '850um')
        self._logger.debug('Done accumulate_bp.')


class JCMTSL450um(JCMTSLSpatialSpectral):

    def accumulate_blueprint(self, bp):
        """Configure the telescope-specific ObsBlueprint at the CAOM model Observation level."""
        self._logger.debug('Begin accumulate_bp.')
        super().accumulate_blueprint(bp)
        # 3" from https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/JCMTSL/scubacat.pdf
        bp.set('Chunk.position.resolution', 3.0)
        # 11" from https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/JCMTSL/scubacat.pdf
        bp.set('Chunk.energy.axis.range.start.val', 444.5)
        bp.set('Chunk.energy.axis.range.end.val', 455.5)
        bp.set('Chunk.energy.bandpassName', '450um')
        self._logger.debug('Done accumulate_bp.')


def mapping_factory(clients, config, reporter, observation, storage_name):
    result = None
    if storage_name.product_type() in [ProductType.AUXILIARY, ProductType.NOISE, ProductType.WEIGHT]:
        result = JCMTSLBase(storage_name, clients, reporter, observation, config)
    else:
        if '850um' in storage_name.file_name:
            result = JCMTSL850um(storage_name, clients, reporter, observation, config)
        else:
            result = JCMTSL450um(storage_name, clients, reporter, observation, config)
    logging.error(f'Created {result.__class__.__name__} for {storage_name.file_uri}')
    return result
