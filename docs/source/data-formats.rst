Data Format Descriptions
========================

TBW
---
TBW consists of a time series of real valued data sampled at f\ :sub:`S` (196 MHz) from
all antennas in the array.  The stand numbering is based on the input into the
digital system rather than the stand number in the array.  The data are divided
into packets that contain either 400 samples (12-bit data) or 1200 samples 
(4-bit data) per stand per polarization.

.. note:: Fields are big endian

+---------+-------------------+---------------------------------+
| Byte(s) | Name              | Description                     |
+=========+===================+=================================+
| 0-3     | Sync Word         | 0xDECODE5C                      |
+---------+-------------------+---------------------------------+
| 4-6     | Frame Count       | Frame number within capture     |
+---------+-------------------+---------------------------------+
| 7       | ID                | Always zero                     |
+---------+-------------------+---------------------------------+
| 8-11    | Second count      | Always zero                     |
+---------+-------------------+---------------------------------+
| 12-13   | TBW ID            | Packed frame source ID          |
+---------+-------------------+---------------------------------+
| 14-15   | Unassigned        | Always zero                     |
+---------+-------------------+---------------------------------+
| 16-23   | Time tag          | Time tag [#F1]_                 |
+---------+-------------------+---------------------------------+
| 24-1223 | Data              | ``i4`` or ``i12`` data - [time,]|
+---------+-------------------+---------------------------------+

The packing on the ``TBW ID`` field is:

+--------+--------------+
| Bit(s) | Description  |
+========+==============+
| 0-13   | Stand        |
+--------+--------------+
| 14     | Bits (0 = 12)|
+--------+--------------+
| 15     | Is TBW?      |
+--------+--------------+

Both the ``i4`` and ``i12`` versions of the data are stored as two's complement
integers.

TBN
---
TBN data is similar to TBW data, except that the data are complex and have a
variable sample rate of up to 100 kHz. The data is divided into packets with 512
samples per stand per polarization.

.. note:: Fields are big endian

+---------+-------------------+---------------------------------+
| Byte(s) | Name              | Description                     |
+=========+===================+=================================+
| 0-3     | Sync Word         | 0xDECODE5C                      |
+---------+-------------------+---------------------------------+
| 4-6     | Frame Count       | Always zero                     |
+---------+-------------------+---------------------------------+
| 7       | ID                | Always zero                     |
+---------+-------------------+---------------------------------+
| 8-11    | Tuning Word       | Tuning frequency [#F2]_         |
+---------+-------------------+---------------------------------+
| 12-13   | TBN ID            | Packed frame source ID          |
+---------+-------------------+---------------------------------+
| 14-15   | Gain              | TBN gain used                   |
+---------+-------------------+---------------------------------+
| 16-23   | Time tag          | Time tag [#F1]_                 |
+---------+-------------------+---------------------------------+
| 24-1047 | Data              | ``ci8`` data - [time,]          |
+---------+-------------------+---------------------------------+

The packing on the ``TBN ID`` field is:

+--------+--------------+
| Bit(s) | Description  |
+========+==============+
| 0-13   | Stand        |
+--------+--------------+
| 14     | Reserved     |
+--------+--------------+
| 15     | Is TBW?      |
+--------+--------------+

``ci8`` is also known as "8+8-bit complex integer".

DRX
---
DRX data consist of a time series of complex data sampled at f\ :sub:`S` / ``Decimation Factor``.
The data are divided into packets of 4096 samples per beam per tuning per
polarization.  Data is typically recorded such that only a single beam is
present.

.. note:: Fields are big endian

+---------+-------------------+---------------------------------+
| Byte(s) | Name              | Description                     |
+=========+===================+=================================+
| 0-3     | Sync Word         | 0xDECODE5C                      |
+---------+-------------------+---------------------------------+
| 4-6     | Frame Count       | Always zero                     |
+---------+-------------------+---------------------------------+
| 7       | ID                | Packed frame source ID          |
+---------+-------------------+---------------------------------+
| 8-11    | Second Count      | Always zero                     |
+---------+-------------------+---------------------------------+
| 12-13   | Decimation Factor | Decimation factor of f\ :sub:`S`|
+---------+-------------------+---------------------------------+
| 14-15   | Time Offset       | Time correction to the time tag |
+---------+-------------------+---------------------------------+
| 16-23   | Time tag          | Time tag [#F1]_                 |
+---------+-------------------+---------------------------------+
| 24-28   | Flags             | Always zero                     |
+---------+-------------------+---------------------------------+
| 29-31   | Tuning Word       | Tuning frequency [#F2]_         |
+---------+-------------------+---------------------------------+
| 32-4127 | Data              | ``ci4`` data - [time,]          |
+---------+-------------------+---------------------------------+

The packing on the ``ID`` field is:

+--------+--------------+
| Bit(s) | Description  |
+========+==============+
| 0-2    | Beam         |
+--------+--------------+
| 3-5    | Tuning       |
+--------+--------------+
| 6      | Reserved     |
+--------+--------------+
| 7      | Polarization |
+--------+--------------+

The ``ci4`` aka "4+4-bit complex integer" data are stored as two's complement
integers with the real part stored in the first four bits.

DR Spectrometer
---------------
DR Spectrometer data contain DRX data that has been transformed to the Fourier
domain, converted to power, and integrated on-the-fly from the data recorder.  
The format allows for variable integration times, channel counts, and
polarization products to be stored.

+---------+-------------------+---------------------------------+
| Byte(s) | Name              | Description                     |
+=========+===================+=================================+
| 0-3     | Magic Word 1      | 0xC0DEC0DE                      |
+---------+-------------------+---------------------------------+
| 4-11    | Time tag          | Time tag of first frame [#F1]_  |
+---------+-------------------+---------------------------------+
| 12-13   | Time offset       | Time correction to the time tag |
+---------+-------------------+---------------------------------+
| 14-15   | Decimation Factor | Decimation factor of f\ :sub:`S`|
+---------+-------------------+---------------------------------+
| 16-23   | Tuning Word[2]    | Tuning frequencies [#F2]_ [#F3]_|
+---------+-------------------+---------------------------------+
| 24-39   | Fills[4]          | Fill counts [#F4]_              |
+---------+-------------------+---------------------------------+
| 40-43   | Errors[4]         | Error flags [#F4]_              |
+---------+-------------------+---------------------------------+
| 44      | Beam              | Beam                            |
+---------+-------------------+---------------------------------+
| 45      | Stokes Format     | Polarization format code        |
+---------+-------------------+---------------------------------+
| 46      | Version           | Packet format version           |
+---------+-------------------+---------------------------------+
| 47      | Flags             | Bit flags                       |
+---------+-------------------+---------------------------------+
| 48-51   | N\ :sub:`freq`    | Transform length [#F5]_         |
+---------+-------------------+---------------------------------+
| 52-53   | N\ :sub:`int`     | Integration count [#F5]_        |
+---------+-------------------+---------------------------------+
| 54-67   | Sats[4]           | Saturation counts [#F4]_        |
+---------+-------------------+---------------------------------+
| 68-71   | Magic Word 2      | 0xED0CED0C                      |
+---------+-------------------+---------------------------------+
| 72-N    | Data              | ``float`` data - [pol,chan]     |
+---------+-------------------+---------------------------------+

The ``Stokes Format`` is a bit field defines what polarizations are included in
the packet data.  The fields are:

+-------+--------------+
| Value | Pol. Product |
+=======+==============+
| 0x01  | XX*          |
+-------+--------------+
| 0x02  | Real(XY*)    |
+-------+--------------+
| 0x04  | Imag(XY*)    |
+-------+--------------+
| 0x08  | YY*          |
+-------+--------------+
| 0x10  | I            |
+-------+--------------+
| 0x20  | Q            |
+-------+--------------+
| 0x40  | U            |
+-------+--------------+
| 0x80  | V            |
+-------+--------------+

The data are stored as little endian ``float`` values.

TBF
---
TBF is similar to both TBW and TBN, but is a complex frequency domain product
that contains blocks of 12 channels from all stands and polarizations.  Each
channel has a bandwidth of f\ :sub:`C` (25 kHz) and there may be up to 132
different values of ``First Channel`` within a single recording.  The stand
ordering is based on the input into the digital system rather than the stand
number in the array.  

.. note:: Fields are big endian

+---------+-------------------+---------------------------------+
| Byte(s) | Name              | Description                     |
+=========+===================+=================================+
| 0-3     | Sync Word         | 0xDECODE5C                      |
+---------+-------------------+---------------------------------+
| 4-6     | Frame Count       | Always zero                     |
+---------+-------------------+---------------------------------+
| 7       | ADP ID            | Always one                      |
+---------+-------------------+---------------------------------+
| 8-11    | Second Count      | Always zero                     |
+---------+-------------------+---------------------------------+
| 12-13   | First Channel     | First channel in packet         |
+---------+-------------------+---------------------------------+
| 14-15   | Unassigned        | Always zero                     |
+---------+-------------------+---------------------------------+
| 16-23   | Time tag          | Time tag [#F1]_                 |
+---------+-------------------+---------------------------------+
| 24-6167 | Data              | ``ci4`` data - [chan,stand,pol] |
+---------+-------------------+---------------------------------+

The ``ci4`` aka "4+4-bit complex integer" data are stored as two's complement
integers with the real part stored in the first four bits.

COR
---
The COR format contains full polarization complex visibility data from a single
baseline pair.  The stand numbering for the baseline pair is based on the input
into the digital system rather than the stand number in the array.  

.. note:: Fields are big endian

+---------+-------------------+-------------------------------------------+
| Byte(s) | Name              | Description                               |
+=========+===================+===========================================+
| 0-3     | Sync Word         | 0xDECODE5C                                |
+---------+-------------------+-------------------------------------------+
| 4-6     | Frame Count       | Always zero                               |
+---------+-------------------+-------------------------------------------+
| 7       | ADP ID            | Always one                                |
+---------+-------------------+-------------------------------------------+
| 8-11    | Second Count      | Always zero                               |
+---------+-------------------+-------------------------------------------+
| 12-13   | First Channel     | First channel in packet                   |
+---------+-------------------+-------------------------------------------+
| 14-15   | Gain              | COR gain used                             |
+---------+-------------------+-------------------------------------------+
| 16-23   | Time tag          | Time tag [#F1]_                           |
+---------+-------------------+-------------------------------------------+
| 24-27   | Navg              | Integration time in samples               |
+---------+-------------------+-------------------------------------------+
| 28-29   | Stand 1           | First stand in baseline                   |
+---------+-------------------+-------------------------------------------+
| 30-31   | Stand 2           | Second stand in baseline                  |
+---------+-------------------+-------------------------------------------+
| 32-4255 | Data              | ``Complex<float>`` data - [chan,pol1,pol2]|
+---------+-------------------+-------------------------------------------+

The data are stored as little endian ``Complex<float>`` values.

.. [#F1] Time tags are expressed as ``uint64_t`` integers in units of ticks of a
         clock at f\ :sub:`S` since the start of the UNIX epoch (1970 Jan 1 00:00:00 UTC)
.. [#F2] Tunings are expressed as ``uint32_t`` integers in units of
         f\ :sub:`S` / 2\ :sup:`32` Hz
.. [#F3] Ordering is [Tuning 1, Tuning 2]
.. [#F4] Ordering is [Tuning 1/pol 0, Tuning 2/pol 0, Tuning 1/pol 1, Tuning 2/pol 1]
.. [#F5] Valid values can be found at https://leo.phys.unm.edu/~lwa/astro/scheds/spec.html
