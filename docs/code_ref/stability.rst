Current status of sub-packages
******************************

SunPy has in variations in stability across sub-packages.
This document summarizes the current status of the SunPy sub-packages, so that users understand where they might expect changes in future, and which sub-packages they can safely use for production code.
For help or to discuss specific sub-packages, refer to the `maintainer list <https://sunpy.org/project/#maintainers>`__ to find who to contact.

The classification is as follows:

{% raw %}
.. raw:: html
{% endraw %}
   <style>
         .planned:before {
              color: #cbcbcb;
              content: "⬤";
              padding: 5px;
         }
         .dev:before {
              color: #ffad00;
              content: "⬤";
              padding: 5px;
         }
         .stable:before {
              color: #4e72c3;
              content: "⬤";
              padding: 5px;
         }
         .mature:before {
              color: #03a913;
              content: "⬤";
              padding: 5px;
         }
         .pendingdep:before {
              color: #a84b03;
              content: "⬤";
              padding: 5px;
         }
         .deprecated:before {
              color: #ff0000;
              content: "⬤";
              padding: 5px;
         }
    </style>

    <table align='center'>
      <tr>
        <td align='center'><span class="planned"></span></td>
        <td>Planned</td>
      </tr>
      <tr>
        <td align='center'><span class="dev"></span></td>
        <td>Actively developed, be prepared for possible significant changes.</td>
      </tr>
      <tr>
        <td align='center'><span class="stable"></span></td>
        <td>Reasonably stable, any significant changes/additions will generally include backwards-compatiblity.</td>
      </tr>
      <tr>
        <td align='center'><span class="mature"></span></td>
        <td>Mature.  Additions/improvements possible, but no major changes planned. </td>
      </tr>
      <tr>
        <td align='center'><span class="pendingdep"></span></td>
        <td>Pending deprecation.  Might be deprecated in a future version.</td>
      </tr>
      <tr>
        <td align='center'><span class="deprecated"></span></td>
        <td>Deprecated.  Might be removed in a future version.</td>
      </tr>
    </table>

The current planned and existing sub-packages are:

{% raw %}
.. raw:: html
{% endraw %}

    <table border="1" class="docutils stability" align='center'>
        <tr>
            <th class="head">
                Sub-Package
            </th>
            <th class="head">
                &nbsp;
            </th>
            <th class="head">
                Comments
            </th>
        </tr>
    {% for module, prop in sunpy_modules.items() %}
        <tr>
            <td>
                <a href="../code_ref/{{ module }}.html">sunpy.{{ module }}</a>
            </td>
            <td align='center'>
                <span class="{{ prop['status'] }}"></span>
            </td>
            <td>
                {{ prop['comments'] }}
            </td>
        </tr>
    {% endfor %}
    </table>


Taken with love from the `Astropy project. <https://github.com/astropy/astropy/blob/master/LICENSE.rst>`_
